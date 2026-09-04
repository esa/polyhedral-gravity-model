# Profiling the Kernels with Nsight Compute

`src/benchmark.cpp` builds into `polyhedralGravity_benchmark`, a driver which evaluates one polyhedron at
`N` random computation points without a Python interpreter in between. It exists so that the kernels can be
handed to a profiler; this file collects the command lines which do that with NVIDIA's Nsight Compute
(`ncu`), so that the resulting `.ncu-rep` can be opened in `ncu-ui`.

Everything below is meant to be copied and pasted. The commands assume the repository root as the working
directory and use the [shell variables](#0-shell-preamble) set up in the preamble.

---

## Prerequisites

* An NVIDIA GPU and a CUDA build of this project (`POLYHEDRAL_GRAVITY_DEVICE_BACKEND=CUDA`). Nsight Compute
  profiles CUDA only; for the HIP backend the equivalent tool is `rocprof`/ `omniperf`.
* `ncu` and `ncu-ui`, both part of the CUDA Toolkit (`/usr/local/cuda/bin`) and of the standalone Nsight
  Compute installation. Check with `ncu --version`.
* **Permission to read the performance counters.** Without it every run fails with
  `ERR_NVGPU_DEBUG_MODE: The user does not have permission to access NVIDIA GPU Performance Counters`.
  Either profile as root (`sudo -E ncu ...`, `-E` keeps the environment), or lift the restriction
  permanently:

  ```bash
  echo 'options nvidia NVreg_RestrictProfilingToAdminUsers=0' | sudo tee /etc/modprobe.d/nvidia-profiler.conf
  sudo update-initramfs -u   # Debian/ Ubuntu; other distributions: rebuild the initramfs accordingly
  sudo reboot
  ```

---

## 1. Build for profiling

`CMakePresets.json` in the repository root carries the configurations this file profiles. Pick the fast math
one, which is what the GPU numbers in the repository are taken with:

```bash
cmake --preset benchmark-cuda-fast-math
cmake --build --preset benchmark-cuda-fast-math
# -> build/benchmark-cuda-fast-math/polyhedralGravity_benchmark
```

| Preset                     | Device backend | Fast math | Use it for                                                              |
|----------------------------|----------------|-----------|--------------------------------------------------------------------------|
| `benchmark-cuda-fast-math` | `CUDA`         | `ON`      | **the default here** — the configuration the recipes below profile        |
| `benchmark-cuda`           | `CUDA`         | `OFF`     | the same kernel under strict IEEE-754, i.e. what fast math actually buys  |
| `benchmark-fast-math`      | `AUTO`         | `ON`      | a machine whose backend is HIP/ SYCL, or a CPU-only host                  |
| `benchmark`                | `AUTO`         | `OFF`     | plain timing runs anywhere, no profiler involved                          |

Every one of them configures `Release`, switches the tests, the Python interface, and the standalone
executable off, and turns two things on:

| Cache variable                        | Why it is in the preset                                                                                       |
|---------------------------------------|-----------------------------------------------------------------------------------------------------------------|
| `BUILD_POLYHEDRAL_GRAVITY_BENCHMARK`  | builds `polyhedralGravity_benchmark` itself (fetches `argparse`)                                                  |
| `POLYHEDRAL_GRAVITY_LINE_INFO`        | adds `-lineinfo`, without which ncu cannot attribute counters to source lines and its Source page stays empty     |
| `CMAKE_POLICY_VERSION_MINIMUM`        | `3.5`, which CMake 4 needs to accept the bundled TetGen's build system. Ignored by older CMake, as in `pyproject.toml` |

`POLYHEDRAL_GRAVITY_FAST_MATH` only ever changes the `FLOAT32` kernel — it becomes `-use_fast_math` for nvcc,
i.e. the hardware approximations for division, square root, and the transcendentals instead of the corrected
IEEE-754 ones. `FLOAT64` is bit-identical with and without it, so a `f64` profile is the same for both
presets. Profiling `benchmark-cuda` against `benchmark-cuda-fast-math` and comparing the instruction mix of
the two reports is the direct way to see where those approximations land.

Two things a preset cannot do for you:

* **`CXX` has to understand CUDA**, since Kokkos compiles the device code with the C++ compiler: a
  CUDA-capable `clang++` or NVIDIA's `nvcc_wrapper` (see the README's build section). Presets do not pin it,
  so pass it in the environment on the first configure of a build directory:

  ```bash
  CXX=/path/to/kokkos/bin/nvcc_wrapper cmake --preset benchmark-cuda-fast-math
  ```

* **The generator** stays whatever your platform defaults to. Add `-G Ninja` on the first configure if you
  want Ninja, or set `CMAKE_GENERATOR` in the environment; both survive into the preset's build directory.

`-lineinfo` costs nothing at runtime, so this build is also the one to take the driver's own timings from.

---

## 2. The meshes to profile with

The driver takes any format the library reads, i.e. a `.node`/`.face` pair (in that order) or a single
`.obj`/`.tab`/`.off`/`.ply`/`.stl`/`.mesh` file. Everything below lives in `examples/data`:

| Mesh                            | Faces   | What it is good for                                                        |
|---------------------------------|---------|----------------------------------------------------------------------------|
| `tsoulis.node tsoulis.face`     | 12      | smoke test only — 12 faces do not fill a single warp                        |
| `cube.node cube.face`           | 12      | same, but with an analytical solution to compare against                    |
| `216Kleopatra.obj`              | 4 092   | a small but realistic body; one team fits its faces in few iterations       |
| `Eros.node Eros.face`           | 14 744  | **the default choice** — the mesh every number in the repository refers to  |
| `Itokawa`, `67P`, `Ryugu`, …    | 200k–3M | the large end; fetch them with `examples/data/download_models.sh`           |

The face count is the inner loop of `evaluateMultiPoint` (one team per point, `TeamThreadRange` over the
faces), the point count `-n` is its grid. Both matter: a mesh with few faces gives short teams, and a small
`-n` leaves the SMs empty regardless of the mesh.

---

## 3. The driver itself

```bash
# a plain timing run, no profiler: the table of best/ mean runtime and the accuracy against a FLOAT64 CPU reference
./build/benchmark-cuda-fast-math/polyhedralGravity_benchmark examples/data/Eros.node examples/data/Eros.face \
    -b gpu cpu -p f32 f64 -n 100000 -r 10
```

Under a profiler the flags to use are different — one launch of the kernel under investigation and nothing
else around it:

| Flag                     | Value under `ncu`        | Why                                                                                          |
|--------------------------|--------------------------|-----------------------------------------------------------------------------------------------|
| `-b, --backend`          | `gpu`                    | `serial`/ `cpu` launch no device kernels at all                                                |
| `-p, --precision`        | `f32` (or `f64`)         | one precision per report unless you deliberately compare, see [snippet 5](#5-f32-against-f64) |
| `-n, --points`           | `100000`                 | enough teams to fill the GPU; the profiler's replay cost grows with it                         |
| `-r, --repetitions`      | `1`                      | every repetition is another launch ncu would replay                                            |
| `--skip-accuracy-check`  | set                      | skips the `CPU_PARALLEL`/`FLOAT64` reference evaluation, which is pure host time               |
| `-s, --seed`             | `42` (default)           | the point sample only depends on the seed and the box, so two runs profile the very same work  |

`--integrity disable` is already the default, i.e. no mesh checking kernels run. Full option list:
`./build/benchmark-cuda-fast-math/polyhedralGravity_benchmark --help`.

Where the points sit is a knob of its own: by default they are drawn from the mesh's bounding box scaled by
`--sample-scale 2.0`, so roughly a tenth of them fall inside the body and the rest around it. Points inside
and points close to a face take the numerically critical branches of `evaluateFace`, which is exactly the
divergence a profile is meant to show. `--sample-scale 1.0` keeps everything within the bounding box,
`--sample-scale 10.0` moves the sample far away from the body, and `--sample-box MIN MAX` replaces the
heuristic altogether.

### What one run launches

Kokkos names its kernels, but ncu identifies them by the C++ symbol, which is the enclosing function of the
Kokkos lambda. The right hand column is what to filter on with `-k regex:`:

| Kokkos label                              | ncu sees (demangled, substring)   | When                                                    |
|-------------------------------------------|-----------------------------------|---------------------------------------------------------|
| `polyhedralGravity::evaluateMultiPoint`   | `runMultiPointKernel`             | once per `evaluable(points)` — **the kernel of interest** |
| `polyhedralGravity::evaluate`             | `runSinglePointKernel`            | only for single point evaluations, not used by the driver |
| `polyhedralGravity::initializeFaceProperties` | `runInitializationKernel`     | once per `GravityEvaluable`, i.e. per backend/ precision pair |
| `polyhedralGravity::narrowVertices`       | `narrowMesh`                      | once per `GravityEvaluable`, `f32` only                  |
| `polyhedralGravity::smallestFaceIndex`    | `shiftFaceIndicesToZeroBased`     | once while the mesh is read                              |

With `-r 1` the driver launches `runMultiPointKernel` **twice**: once untimed to warm up, once timed. That
is what `--launch-skip 1 --launch-count 1` below refers to — profile the second, steady-state launch.

---

## 4. Copy + paste

### 0. Shell preamble

```bash
# from the repository root
export PG_BENCH="$PWD/build/benchmark-cuda-fast-math/polyhedralGravity_benchmark"
export PG_DATA="$PWD/examples/data"
export PG_REPORTS="$PWD/profiling-results" && mkdir -p "$PG_REPORTS"   # .gitignore keeps the reports out of git
export CUDA_VISIBLE_DEVICES=0                            # pick one GPU on a multi GPU machine

# An array, so that a two file mesh stays two arguments. A single file mesh is a single entry:
# PG_MESH=("$PG_DATA/216Kleopatra.obj")
PG_MESH=("$PG_DATA/Eros.node" "$PG_DATA/Eros.face")
```

### 1. Which kernels and how many launches?

Cheap first pass. It profiles everything, prints one line per kernel, and writes no report — run it once so
that the launch indices in the snippets below are the ones you expect.

```bash
ncu --section LaunchStats --print-summary per-kernel \
    "$PG_BENCH" "${PG_MESH[@]}" -b gpu -p f32 -n 100000 -r 1 --skip-accuracy-check
```

### 2. The report to open in `ncu-ui`

The one to reach for by default: every section on the timed launch of the multi point kernel, with the
sources embedded so that the report can be opened on a machine which does not have the checkout.

```bash
ncu --set full \
    --kernel-name-base demangled -k regex:runMultiPointKernel \
    --launch-skip 1 --launch-count 1 \
    --import-source yes \
    --force-overwrite --export "$PG_REPORTS/eros_f32_full" \
    "$PG_BENCH" "${PG_MESH[@]}" -b gpu -p f32 -n 100000 -r 1 --skip-accuracy-check

ncu-ui "$PG_REPORTS/eros_f32_full.ncu-rep"
```

`--set full` collects everything including the roofline chart and replays the kernel a few dozen times; on
100 000 points expect a couple of minutes. Drop to `--set detailed` (no source counters) or
`--set default` if that is too slow.

### 3. Roofline only

```bash
ncu --section SpeedOfLight --section SpeedOfLight_RooflineChart \
    --kernel-name-base demangled -k regex:runMultiPointKernel \
    --launch-skip 1 --launch-count 1 \
    --force-overwrite --export "$PG_REPORTS/eros_f32_roofline" \
    "$PG_BENCH" "${PG_MESH[@]}" -b gpu -p f32 -n 100000 -r 1 --skip-accuracy-check
```

Open it and go to *Details → GPU Speed Of Light Throughput → Roofline*. The kernel reads a handful of
vertices per face and then does dozens of divisions, square roots, and logarithms on them, so it should sit
far to the right of the ridge point, i.e. compute/ latency bound rather than memory bound. If a change moves
it left, the change made the arithmetic cheaper or the memory traffic worse.

### 4. Source level hot spots

Needs the `POLYHEDRAL_GRAVITY_LINE_INFO=ON` which every benchmark preset sets.

```bash
ncu --section SourceCounters --section Occupancy --section LaunchStats \
    --kernel-name-base demangled -k regex:runMultiPointKernel \
    --launch-skip 1 --launch-count 1 \
    --import-source yes \
    --force-overwrite --export "$PG_REPORTS/eros_f32_source" \
    "$PG_BENCH" "${PG_MESH[@]}" -b gpu -p f32 -n 100000 -r 1 --skip-accuracy-check
```

In `ncu-ui` this is the *Source* page: instructions executed and stall reasons per line of
`GravityModelDetail.h`/ `GravityEvaluable.cpp`, which is where `evaluateFace` shows what its branches cost.

### 5. `f32` against `f64`

Both precisions in a single report. The driver builds one evaluable per precision, so the four matching
launches are, in order: `f32` warm-up, `f32` timed, `f64` warm-up, `f64` timed. `ncu-ui` lists them
separately and lets you set any of them as the baseline for the others.

```bash
ncu --set full \
    --kernel-name-base demangled -k regex:runMultiPointKernel \
    --launch-count 4 \
    --import-source yes \
    --force-overwrite --export "$PG_REPORTS/eros_f32_vs_f64" \
    "$PG_BENCH" "${PG_MESH[@]}" -b gpu -p f32 f64 -n 100000 -r 1 --skip-accuracy-check
```

### 6. The occupancy numbers, printed to the terminal

The counters behind the `Kokkos::LaunchBounds<128, 6>` in `runMultiPointKernel`: registers per thread,
achieved occupancy, and the local memory traffic which is what a register spill looks like. No report file,
just a table — the fastest way to check whether a change to the kernel body pushed it over the cliff.

```bash
ncu --metrics launch__registers_per_thread,\
launch__occupancy_limit_registers,\
sm__warps_active.avg.pct_of_peak_sustained_active,\
smsp__sass_inst_executed_op_local_ld.sum,\
smsp__sass_inst_executed_op_local_st.sum,\
gpu__time_duration.sum \
    --kernel-name-base demangled -k regex:runMultiPointKernel \
    --launch-skip 1 --launch-count 1 \
    "$PG_BENCH" "${PG_MESH[@]}" -b gpu -p f32 -n 100000 -r 1 --skip-accuracy-check
```

Non-zero `local_ld`/`local_st` means the kernel spills; that is the regression the launch bounds comment in
`GravityEvaluable.cpp` warns about. Add `--csv` to paste the output into a spreadsheet, and use
`ncu --query-metrics` to check what else this architecture exposes.

### 7. A sweep over meshes and precisions

```bash
for MESH_NAME in 216Kleopatra Eros; do
  case "$MESH_NAME" in
    Eros) MESH=("$PG_DATA/Eros.node" "$PG_DATA/Eros.face") ;;
    *)    MESH=("$PG_DATA/$MESH_NAME.obj") ;;
  esac
  for PRECISION in f32 f64; do
    ncu --set full \
        --kernel-name-base demangled -k regex:runMultiPointKernel \
        --launch-skip 1 --launch-count 1 \
        --import-source yes \
        --force-overwrite --export "$PG_REPORTS/${MESH_NAME}_${PRECISION}" \
        "$PG_BENCH" "${MESH[@]}" -b gpu -p "$PRECISION" -n 100000 -r 1 --skip-accuracy-check
  done
done
```

### 8. Before and after a change to the kernel

Profile the old code, keep the report, profile the new one, and compare the two in the UI:

```bash
git stash                                   # or: git checkout <the commit to compare against>
cmake --build --preset benchmark-cuda-fast-math -j "$(nproc)"
ncu --set full --kernel-name-base demangled -k regex:runMultiPointKernel \
    --launch-skip 1 --launch-count 1 --import-source yes \
    --force-overwrite --export "$PG_REPORTS/baseline" \
    "$PG_BENCH" "${PG_MESH[@]}" -b gpu -p f32 -n 100000 -r 1 --skip-accuracy-check

git stash pop
cmake --build --preset benchmark-cuda-fast-math -j "$(nproc)"
ncu --set full --kernel-name-base demangled -k regex:runMultiPointKernel \
    --launch-skip 1 --launch-count 1 --import-source yes \
    --force-overwrite --export "$PG_REPORTS/candidate" \
    "$PG_BENCH" "${PG_MESH[@]}" -b gpu -p f32 -n 100000 -r 1 --skip-accuracy-check

# open the old one first, press "Add Baseline", then open the new one on top of it
ncu-ui "$PG_REPORTS/baseline.ncu-rep" "$PG_REPORTS/candidate.ncu-rep"
```

Every section of the candidate then shows the relative change against the baseline next to each metric. Read a runtime *number* off the driver's
own table instead (`-r 10`, no profiler): under ncu the kernel is replayed and its clocks are fixed, so the
durations in the report compare well against each other but not against a normal run.

### 9. Fast math against IEEE-754

What `-use_fast_math` does to the `FLOAT32` kernel, in one comparison: configure both presets once, then
profile the same launch out of each build directory.

```bash
cmake --preset benchmark-cuda            && cmake --build --preset benchmark-cuda            -j "$(nproc)"
cmake --preset benchmark-cuda-fast-math  && cmake --build --preset benchmark-cuda-fast-math  -j "$(nproc)"

for BUILD in benchmark-cuda benchmark-cuda-fast-math; do
  ncu --set full \
      --kernel-name-base demangled -k regex:runMultiPointKernel \
      --launch-skip 1 --launch-count 1 \
      --import-source yes \
      --force-overwrite --export "$PG_REPORTS/$BUILD" \
      "$PWD/build/$BUILD/polyhedralGravity_benchmark" "${PG_MESH[@]}" \
      -b gpu -p f32 -n 100000 -r 1 --skip-accuracy-check
done

ncu-ui "$PG_REPORTS/benchmark-cuda.ncu-rep" "$PG_REPORTS/benchmark-cuda-fast-math.ncu-rep"
```

The difference to look for is in *Instruction Statistics* and on the Source page: the corrected division and
square root sequences collapse into the hardware approximations, which is where the roughly 30% the option
is documented to buy come from. Run the same two binaries without `ncu` (`-r 10`) for the wall clock number,
and note that the accuracy column of the driver's own table is the other half of that trade.

---

## Reading a report without the UI

```bash
ncu --import "$PG_REPORTS/eros_f32_full.ncu-rep" --page details        # what the Details page shows
ncu --import "$PG_REPORTS/eros_f32_full.ncu-rep" --page details --csv  # the same, machine readable
ncu --import "$PG_REPORTS/eros_f32_full.ncu-rep" --page raw            # every collected metric
```

---

## Pitfalls

* **Profiling the warm-up instead of the real launch.** The first `evaluateMultiPoint` of a run pays for the
  lazily created kernel and a cold cache. `--launch-skip 1 --launch-count 1` is what skips it.
* **Forgetting `-k`.** Without a kernel filter ncu replays the initialization and mesh kernels as well; with
  `--set full` that turns a two minute run into a long one for data nobody looks at.
* **Timings taken under the profiler.** `ncu` fixes the clocks (`--clock-control base`) and replays each
  kernel several times. Compare reports with reports and wall clock times with wall clock times. Pass
  `--clock-control none` if you specifically want the kernel timed at the GPU's own boost behaviour.
* **`-n` too small.** With a few thousand points the grid does not fill the GPU and every occupancy number
  describes the tail, not the kernel. 100 000 points on the Eros mesh is a sensible working point.
* **An empty Source page.** The build was configured without `POLYHEDRAL_GRAVITY_LINE_INFO=ON` (all four
  presets set it, a hand written `cmake ..` does not), or the
  report was moved to another machine without `--import-source yes`.
* **`--set full` on a 3M face mesh.** The replay cost scales with the kernel, so start from `--set default`
  or reduce `-n` before reaching for the large bodies.

---

## Related: a timeline instead of a kernel

`ncu` answers "why is this kernel slow"; `nsys` answers "what does the whole run spend its time on",
including the host side, the point upload, and the result download:

```bash
nsys profile --trace=cuda,nvtx,osrt \
     --force-overwrite true --output "$PG_REPORTS/eros_timeline" \
     "$PG_BENCH" "${PG_MESH[@]}" -b gpu -p f32 -n 100000 -r 10 --skip-accuracy-check

nsys-ui "$PG_REPORTS/eros_timeline.nsys-rep"
```

Note the different flags: `-r 10` and no launch filtering, since a timeline is about the sequence of
launches and copies rather than about one of them.
