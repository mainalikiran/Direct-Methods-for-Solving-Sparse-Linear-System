
function xsuper = supernode(Ap,parent)
% p = amd(A);
% Ap = A(p,p);
% parent = etree(Ap);
nsuper = 0;
n = size(Ap,1);
sizes = treesize(parent);
children = childrencount(parent);    
xsuper = [];
prevnonz = zeros(1,n);
for j = 1:n
    isleaf=0;
    idx = find(Ap(j+1:n,j)) + j;
    for p = 1:size(idx,1)
        t = idx(p);
        k = prevnonz(t);
        if k < j- sizes(j) +1
            isleaf=1;
        end
        prevnonz(t) = j;
    end
    if isleaf == 1 || children(j) > 1
        nsuper = nsuper+1;
        xsuper(nsuper) = j;
    end  
end
xsuper(nsuper+1)=n+1;
end
