function [supidx] = toSupidx( columns, supmembership,xsuper )

scols = size(columns);
supidx = [];

if min(scols)>0

marker = zeros(max(size(xsuper)),1);

dim = 1;
if scols(2) > 1
    dim =2;
end

q = 1;
while q<=scols(dim)
    c = columns(q);
    super =  supmembership(c);
    dimk = xsuper(super+1) - xsuper(super);
    
    if marker(super)== 0
        supidx = [supidx; super];
        marker(super) = 1;
    end
    
    q = q +1;%dimk;
end

end
end