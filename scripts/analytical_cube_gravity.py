from sympy import *

"""
This python script computes the analytical solution for the potential and acceleration for an arbitrary point P(X,Y,Z)
around _a cube centered at the origin with edge length 2.
"""

init_printing()

# Constants
GRAVITATIONAL_CONSTANT = 6.67430e-11
DENSITY = 1.0

# The point to evaluate, replace values here
X = 1
Y = 1
Z = 0

# Some definitions
_a, _b, _c = symbols('_a _b c')
_x1, _x2, _x3 = symbols('x y z')
_v = _x1 * _x2 * _x3
_r = sqrt(_x1 * _x1 + _x2 * _x2 + _x3 * _x3)
_INNER_SUM_POTENTIAL = ((_v / _a) * log(_a + _r) - (_a * _a / 2) * atan(_v / (_a * _a * _r)))
_INNER_SUM_ACCELERATION = _b * log(_c + _r) + _c * log(_b + _r) - _a * atan(_v / (_a * _a * _r))


def evaluate_potential(X: float, Y: float, Z: float):
    """
    Evaluates the potential at _a given point P(X, Y, Z) for _a cube centered at the origin and edge length 2.
    Args:
        X: first coordinate
        Y: second coordinate
        Z: third coordinate

    Returns:
        potential

    """
    total_sum = _INNER_SUM_POTENTIAL.subs({_a: _x1}) + _INNER_SUM_POTENTIAL.subs({_a: _x2}) + _INNER_SUM_POTENTIAL.subs(
        {_a: _x3})
    total_sum = total_sum.evalf(subs={_x1: 1 - X}) - total_sum.evalf(subs={_x1: -1 - X})
    total_sum = total_sum.evalf(subs={_x2: 1 - Y}) - total_sum.evalf(subs={_x2: -1 - Y})

    # The limit formulation handles the case: 0 * log(..) resulting in zero better than the subs functionality
    if 1 - Z == 0 or -1 - Z == 0:
        total_sum_1 = total_sum.limit(_x3, 1 - Z)
        total_sum_2 = total_sum.limit(_x3, -1 - Z)
        total_sum = (total_sum_1 - total_sum_2).evalf()
    else:
        total_sum = total_sum.evalf(subs={_x3: 1 - Z}) - total_sum.evalf(subs={_x3: -1 - Z})
    return total_sum * GRAVITATIONAL_CONSTANT * DENSITY


def _evaluate_acceleration(X: float, Y: float, Z: float, i: int):
    """
    Evaluates the acceleration at _a given point P(X, Y, Z) for _a cube centered at the origin and edge length 2.
    Args:
        X: first coordinate
        Y: second coordinate
        Z: third coordinate
        i: the cartesian direction

    Returns:
        acceleration for direction i

    """
    total_sum = None
    if i == 1:
        total_sum = _INNER_SUM_ACCELERATION.subs({_a: _x1, _b: _x2, _c: _x3})
    elif i == 2:
        total_sum = _INNER_SUM_ACCELERATION.subs({_a: _x2, _b: _x1, _c: _x3})
    elif i == 3:
        total_sum = _INNER_SUM_ACCELERATION.subs({_a: _x3, _b: _x2, _c: _x1})
    total_sum = total_sum.evalf(subs={_x1: 1 - X}) - total_sum.evalf(subs={_x1: -1 - X})
    total_sum = total_sum.evalf(subs={_x2: 1 - Y}) - total_sum.evalf(subs={_x2: -1 - Y})

    # The limit formulation handles the case: 0 * log(..) resulting in zero better than the subs functionality
    if 1 - Z == 0 or -1 - Z == 0:
        total_sum_1 = total_sum.limit(_x3, 1 - Z)
        total_sum_2 = total_sum.limit(_x3, -1 - Z)
        total_sum = (total_sum_1 - total_sum_2).evalf()
    else:
        total_sum = total_sum.evalf(subs={_x3: 1 - Z}) - total_sum.evalf(subs={_x3: -1 - Z})
    return total_sum * GRAVITATIONAL_CONSTANT * DENSITY * -1.0


def evaluate_acceleration(X: float, Y: float, Z: float):
    """
    Evaluates the acceleration at _a given point P(X, Y, Z) for _a cube centered at the origin and edge length 2.
    Args:
        X: first coordinate
        Y: second coordinate
        Z: third coordinate
        i: the cartesian direction

    Returns:
        acceleration in the three cartesian directions

    """
    return _evaluate_acceleration(X, Y, Z, 1), _evaluate_acceleration(X, Y, Z, 2), _evaluate_acceleration(X, Y, Z, 3)


if __name__ == '__main__':
    print(evaluate_potential(X, Y, Z))
    print(evaluate_acceleration(X, Y, Z))
