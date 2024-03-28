import analytical_cube_gravity
import multiprocessing as mp
import numpy as np
from tqdm import tqdm


def worker(args):
    x, y, z = args
    potential = analytical_cube_gravity.evaluate_potential(x, y, z)
    acceleration = analytical_cube_gravity.evaluate_acceleration(x, y, z)
    return args, potential, acceleration

def main():
    X = np.arange(-5, 5.5, 0.5)
    Y = np.arange(-5, 5.5, 0.5)
    Z = np.arange(-5, 5.5, 0.5)

    computation_points = np.array(np.meshgrid(X, Y, Z)).T.reshape(-1, 3)

    with mp.Pool(mp.cpu_count()) as pool:
        results = list(tqdm(pool.imap(worker, computation_points), total=len(computation_points)))
    with open("../test/resources/analytic_cube_solution_density1.txt", "w") as file:
        for point, potential, acceleration in results:
            print("{:.2f} {:.2f} {:.2f} {} {} {} {}".format(point[0], point[1], point[2], potential, acceleration[0],
                                                            acceleration[1], acceleration[2]), file=file)


if __name__ == '__main__':
    main()
