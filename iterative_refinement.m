%Iterative refinement
%Once we compress the block in ccholblock.m coefficient matrix is not exact
% as it was in the begining so we improve the solution by iterative
% refinement. Here, "A" is the randspd matrix , "b" is block size, "T" is
% tolerance for myrank function, smaller the tolerence lesser the
% compression and "tol" is tolerance for iterative refinement.

function [x ,error] = iterative_refinement(A,b,T,tol)                               
[A2,L,~] = ccholblock(A,b,T) ; %A2 is the input matrix.
n = size(L,1);
x_true = rand(n,1);
c = A2*x_true;
x = L'\(L\c);
r = c - A2 * x;
% r_0 = r;
while norm(r,'fro') >= tol
    d = A2\r;                     %Improve factor
    x = x + d;
    r = c - A2 * x; 
end
error = norm(x-x_true)/norm(x);
% semilogy(r_0,'r');
% legend('residue @begining');
% hold on
% plot(r, 'm');
% semilogy('improved residue')
% hold off
end


