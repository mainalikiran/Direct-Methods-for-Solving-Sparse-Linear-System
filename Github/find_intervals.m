function [intervals] = find_intervals(intervals,frk,lrk)
found = 0;
rowblk=1;
nbd0 = frk;
nbd1 = lrk;
while rowblk<=size(intervals,2)
    
    blkdesc = intervals(:,rowblk);
    bd0 = blkdesc(1);
    bd1 = blkdesc(2);
    if nbd1 < bd0-1
        break;
    end
    if ~(nbd1 < bd0-1 || nbd0 > bd1+1)
        nbd0 = min(nbd0,bd0);  nbd1 = max(nbd1,bd1);
        if found==0
            intervals(:,rowblk) = [nbd0;nbd1];
        else
            %interval needs to be removed
            assert(rowblk>1);
            intervals(:,rowblk-1) = [nbd0;nbd1];
            before = intervals(:,1:rowblk-1);
            if rowblk+1 <= size(intervals,2)
                after = intervals(:,rowblk+1:size(intervals,2));
            else
                after = [];
            end
            intervals = [before,after];
            rowblk=rowblk-1;
        end
        found = 1;
    end
    rowblk=rowblk+1;
    
end
                %block hasnt been found
                if found == 0
                   %find where to insert
                   rowblk=1;
                   while rowblk<=size(intervals,2)
                       if intervals(1,rowblk)>lrk
                           found = 1;
                           before = intervals(:,1:rowblk-1);
                           after = intervals(:,rowblk:size(intervals,2));
                           intervals = [before,[frk ; lrk],after];
                           break;
                       end
                       rowblk=rowblk+1;
                   end
                   if found == 0
                    intervals = [intervals ,[frk ; lrk]];
                   end
                end
end