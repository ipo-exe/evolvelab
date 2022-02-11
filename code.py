import numpy as np


def decode_binary_float(gene, lo=0, hi=1):
    """
    Decode a binary gene
    :param gene: numpy array of binary data (1 and 0)
    :param lo: float of lower boundary
    :param hi: float of upper boundary
    :return: decoded float
    """
    _ids = np.arange(start=0, stop=len(gene))
    _twos = (_ids * 0) + 2
    return ((np.sum(gene * np.power(_twos, _ids)) / np.power(2, len(gene))) * (hi - lo)) + lo


