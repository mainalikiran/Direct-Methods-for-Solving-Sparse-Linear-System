
function [blocksize ,blockindex] = rowstruct(idx)
blocksize = [];
blockindex =[];
b = 1;
index =1;
for j = 1:size(idx,1)
   if j==size(idx,1)
      m1=b;
      blocksize = [blocksize ; b];
      index = idx(j)+1-b;
      blockindex = [blockindex ; index];
  
   elseif (idx(j+1) - idx(j)) == 1
     b = b+1;
    
     continue  
     
     elseif idx(j+1)-idx(j) >1 && (j>=size(idx,1)-1||idx(j+2)-idx(j+1) >=1)
     m2 = b ;
     blocksize = [blocksize ; b];
     index = idx(j)+1-b;
     blockindex = [blockindex ; index];
     b=1;

  else
      b = 1 ;
      m3 = b;
      blocksize = [blocksize ; b];
      index = idx(j);
      blockindex = [blockindex ; index];
   end
end
end