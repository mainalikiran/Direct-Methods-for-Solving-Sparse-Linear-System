%Supernode membership function. Which needs output  "xsuper"  from
%supernode.m and size 'n' of matrix Ap.
%Once We detect the supernode of Ap, this function finds the column index
%lying in corresponding supernode.
% xsuper is size of nsuper. xsuper(I) is first column of supernode I.
% Supmembership is size of n and supmembership(j) is supernode index of
% column j.


function supmembership = member(xsuper ,n)
supmembership = zeros(1,n);
for s = 1:size(xsuper,2)-1
    fc = xsuper(s);
    nc = xsuper(s+1)- xsuper(s);
    for c = fc:fc+nc-1
        supmembership(c) = s;
    end
end