#!/bin/python3

import numpy as np
import numpy.linalg as nla
import scipy.sparse as sp
import scipy.sparse.linalg as sla


def read_mat(fname: str):
    IJV = np.fromfile(fname, sep=" ").reshape(-1, 3)
    I = IJV[:, 0].astype(int) - 1
    J = IJV[:, 1].astype(int) - 1
    return sp.coo_matrix((IJV[:, 2], (I, J))).tocsr()


def read_vec(fname: str):
    return np.fromfile(fname, sep="\n")


def fw_solve(L, b):
    return sla.spsolve(L, b)


if __name__ == "__main__":
    L = read_mat("L.txt")
    U = read_mat("U.txt")
    b = read_vec("b.txt")
    y = read_vec("y.txt")
    ye = fw_solve(L, b)
    err = nla.norm(ye - y)
    print(f"ILU forward solve error: {err}")
    assert nla.norm(ye - y) < 1e-8
