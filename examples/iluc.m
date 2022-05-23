load 'L.txt'
load 'U.txt'
load 'LL.txt'
load 'UU.txt'

L = spconvert(L);
U = spconvert(U);
n = size(L, 1);
I = speye(n);
A = L + U - I;

LL = spconvert(LL);
UU = spconvert(UU);
err = norm(A - LL * UU, Inf);
printf('LU error: %f\n', err);
assert(err < 1e-8);
