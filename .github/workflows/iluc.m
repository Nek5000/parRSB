load 'A.txt'
load 'L.txt'
load 'U.txt'

A = spconvert(A);
L = spconvert(L);
U = spconvert(U);
err = norm(A - L * U, Inf);
printf('LU error: %f\n', err);
assert(err < 1e-8);
