load 'A.txt';
load 'L.txt';
load 'U.txt';

A = spconvert(A);
L = spconvert(L);
U = spconvert(U);

[Lo, Uo] = ilu(A);

Lerr = norm(Lo - L, Inf);
Uerr = norm(Uo - U, Inf);
printf('L err = %f, U err = %f', Lerr , Uerr);
assert(Lerr < 1e-8 && Uerr < 1e-8);
