% finding the vector 's' that counts size of subtree(Number of children) 
%rooted at nodes(vertex) of the graph of the matrix.

function s = treesize(parent)
n = size(parent,2);
s  = ones(1,n);
for i = 1:n-1
    parent1 = parent(i);
    s(parent1) = s(parent1)+ s(i);
end
end