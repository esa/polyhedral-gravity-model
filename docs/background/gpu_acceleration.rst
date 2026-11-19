.. _gpu-acceleration:

GPU Acceleration in v4.0
========================

Version :code:`v4.0` replaces the former TBB/ OpenMP parallelization with a single implementation
written in `Kokkos <https://kokkos.org/>`__. The very same kernel now runs on every
:code:`ComputeBackend` the library was compiled with: :code:`CPU_SERIAL` (Kokkos :code:`Serial`),
:code:`CPU_PARALLEL` (Kokkos :code:`OpenMP`), and - new in :code:`v4.0` - :code:`GPU_PARALLEL`
(Kokkos :code:`Cuda`, :code:`HIP`, or :code:`SYCL`). Independently of the backend, the evaluation can
be carried out in single (:code:`FLOAT32`) or double (:code:`FLOAT64`) precision.

This section presents the resulting runtime on two machines and derives a recommendation on how to
call the model. For how to select a backend and a precision, refer to the
**INSTALLATION & QUICK START** section.

Benchmark Setup
---------------

The measurements were taken with :code:`scripts/measure_runtime.py` and the Python interface
:code:`v4.0.0rc1`, built with the :code:`Serial`, :code:`OpenMP`, and :code:`Cuda` backends and
:code:`POLYHEDRAL_GRAVITY_FAST_MATH=ON` (which only affects :code:`FLOAT32`).
Each run evaluates :math:`N = 1000` random computation points, drawn uniformly from
:math:`[-2, 2]^3`, for the mesh of Eros (24235 vertices and 14744 faces).
The best of three repeats is reported as runtime per point (the lower the better/ faster).

For every backend and both precisions, three ways of calling the model are compared:

- :code:`evaluate` :math:`N \times 1` - the free function, called once per point
- :code:`GravityEvaluable` :math:`N \times 1` - one evaluable, called once per point
- :code:`GravityEvaluable` :math:`1 \times N` - one evaluable, called once with all points as :math:`(N, 3)`-array

The pale bars mark the discouraged ways of calling the model, the saturated bars the recommended one.
Blue denotes :code:`FLOAT32`, red :code:`FLOAT64`.

Workstation: Intel i7-14700K & NVIDIA RTX 5080
----------------------------------------------

.. image:: /figures/runtime_measurements_i7-14700K_RTX5080.png
  :align: center
  :width: 100%
  :alt: Runtime per point for the mesh of (433) Eros on an Intel i7-14700K and an NVIDIA RTX 5080 in v4.0

Superchip: NVIDIA GH200 (Grace CPU & Hopper GPU)
------------------------------------------------

.. image:: /figures/runtime_measurements_GH200.png
  :align: center
  :width: 100%
  :alt: Runtime per point for the mesh of (433) Eros on an NVIDIA GH200 in v4.0

Takeaways
---------

**The GPU pays off for batches of points.** Handing all points over at once to a
:code:`GravityEvaluable` on :code:`GPU_PARALLEL` is by far the fastest configuration on both machines:

.. list-table:: :code:`GravityEvaluable` :math:`1 \times N`, runtime per point in :math:`\mu s`
   :header-rows: 1

   * - Machine
     - Precision
     - :code:`CPU_SERIAL`
     - :code:`CPU_PARALLEL`
     - :code:`GPU_PARALLEL`
   * - i7-14700K & RTX 5080
     - :code:`FLOAT32`
     - 1110.2
     - 72.5
     - 0.7
   * - i7-14700K & RTX 5080
     - :code:`FLOAT64`
     - 974.0
     - 76.9
     - 22.3
   * - GH200
     - :code:`FLOAT32`
     - 2268.4
     - 157.5
     - 0.9
   * - GH200
     - :code:`FLOAT64`
     - 1454.1
     - 95.9
     - 1.7

Compared to :code:`CPU_PARALLEL`, this is a speedup of roughly :math:`100\times` (RTX 5080) and
:math:`175\times` (GH200) in :code:`FLOAT32`, and of :math:`3\times` (RTX 5080) and :math:`55\times` (GH200)
in :code:`FLOAT64`.

**Double precision depends on the GPU.** Consumer GPUs like the RTX 5080 only offer a fraction of their
single precision throughput in double precision, so :code:`FLOAT64` is about :math:`30\times` slower than
:code:`FLOAT32` there. Data center GPUs like the Hopper GPU of the GH200 have dedicated double precision
units and lose less than a factor of two. Keep in mind that :code:`FLOAT32` comes with a relative accuracy
of only about :math:`10^{-4}`, since *Tsoulis et al.*'s formulation cancels large terms against each other.
Prefer :code:`FLOAT64` whenever the result matters more than the throughput.

**Calling the GPU point by point is expensive.** Every call transfers data to the device and launches a
kernel. With a :code:`GravityEvaluable`, whose mesh and cached properties stay on the device, this overhead
amounts to roughly :math:`10` to :math:`50\, \mu s` per point. The free function :code:`evaluate` additionally
uploads and prepares the mesh on every call, which makes it the slowest option on the GPU.
On the many-core Grace CPU, the same per-call overhead even makes :code:`CPU_PARALLEL` slower than
:code:`CPU_SERIAL` for single points.

**Recommendation.** As in :code:`v3.0` (see :ref:`evaluable-vs-eval`), use a :code:`GravityEvaluable`.
Whenever possible, collect the points and evaluate them in one call - especially on the GPU.
For single points, e.g. in a trajectory propagation, :code:`CPU_SERIAL` or :code:`CPU_PARALLEL`
are competitive or even faster, depending on the machine.
