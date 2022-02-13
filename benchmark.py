"""
Module to store benchmark functions

"""
import numpy as np


def griewank_2d(x, y, x0=50, y0=50, level=4):
    """
    compute the 2D Griewank Function (upside down by level)
    :param x: float of x
    :param y: float of y
    :param x0: float of x center
    :param y0: float of y center
    :param level: level to upside down
    :return: float of Griewank Function
    """
    # set center
    x = x - x0
    y = y - y0
    return level - (((np.square(x) + np.square(y)) / 4000) - (np.cos(x) * np.cos(y / np.sqrt(2))) + 1)


def sphere_2d(x, y, x0=50, y0=50, level=100):
    """
    (upside down by level)
    :param x:
    :param y:
    :param x0:
    :param y0:
    :param level:
    :return:
    """
    # set center
    x = x - x0
    y = y - y0
    return level - (np.square(x) + np.square(y))


def rastrigin_2d(x, y, x0=50, y0=50, level=100):
    # set center
    x = x - x0
    y = y - y0
    return level - (20 +
                    (np.square(x) - 10 * np.cos(2 * np.pi * x))
                    + (np.square(y) - 10 * np.cos(2 * np.pi * y)))

def himmelblaus(x, y, x0=50, y0=50, level=1000):
    # set center
    x = x - x0
    y = y - y0
    return level - (np.square(np.square(x) + y - 11) + np.square(x + np.square(y) - 7))
