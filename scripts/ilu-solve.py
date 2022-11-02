#!/bin/python3

import argparse

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


def solve(L, b):
    return sla.spsolve(L, b)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Test ILU solve")
    parser.add_argument(
        "--test_fw", type=int, nargs="?", help="Test forward solve"
    )
    parser.add_argument(
        "--test_bw", type=int, nargs="?", help="Test backward solve"
    )
    args = parser.parse_args()

    L = read_mat("L.txt")
    U = read_mat("U.txt")
    b = read_vec("b.txt")
    y = read_vec("y.txt")
    x = read_vec("x.txt")

    ye = solve(L, b)
    if args.test_fw:
        errf = nla.norm(ye - y)
        print(f"ILU forward solve error: {errf}")
        assert errf < 1e-8

    if args.test_bw:
        r = b - L * U * x
        errb = nla.norm(r[:-1])
        print(f"ILU error: {errb}")
        assert errb < 1e-8
