load 'pre.txt';
load 'post.txt';
A = spconvert(pre);
B = spconvert(post);

[L, U] = ilu(A);
n = size(A, 1);
I = speye(n);
err = norm(L + U - B - I, Inf);
printf('ILU err = %f', err);
