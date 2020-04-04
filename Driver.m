% Finale Driver for result comparision. It compares the results between
% left looking Cholesky decomposition with modified developed algorithm for
% rank sorting Block Low Rank Decomposition applied to left looking
% Cholesky.

% Written By : Mainali, Kiran Kumar
% E-mail: kirankumar.mainali@mavs.uta.edu
% Date: 05/2018
% This is Internship project at Lawrence Berkeley National Laboratory,
% California in association with research scientist at Computational
% Research Devision (CRD), Dr. Mathias Jacquelin.

load aft01
A = Problem.A;
n = size(A,1);
%p = dissect(A);
p = amd(A);

fullrank=0;
T = 1e-6;
[~,~,parent,post,L_symbol] = symbfact(A(p,p),'lo','lower');
p = p(post);
%etreeplot(A(p,p));

parent2 = zeros(1,n);
invpost(post) = 1:n;
for i=1:n
    orig = post(i);
    orig_parent = parent(orig);
    if orig_parent~=0
        newparent = invpost(orig_parent);
    else
        newparent=0;
    end
    parent2(i) = newparent;
end
parent = parent2;

L_symbol = L_symbol(post,post);

%[~,~,parent,post,L_symbol] = symbfact(A(p,p),'lo','lower');
[p2,post2] = etree(A(p,p));
assert(isequal(post2,1:n));
assert(isequal(parent2,p2));
%assert(isequal(parent,p2));

parent = fixtree(parent);
%figure; treeplot(parent);

Ap = A(p,p);

xtrue = rand(n,1);
rhs = Ap*xtrue;

xsuper = supernode(Ap, parent);
%xsuper = supernode(L_symbol, parent);
supmembership = member(xsuper,n);
nsuper = size(xsuper,2);
L_exact = chol(Ap,'lower');

L1 = tril(Ap);
results1 = containers.Map('KeyType','char','ValueType','any');
%this is a full rank factorization written from scratch to double check
%VERY SLOW leave if 0 to not use it

if false
    L2 = tril(Ap);
    for i = 1:nsuper-1
        dimi = xsuper(i+1)- xsuper(i);
        fc = xsuper(i);
        
        for tgt=fc:fc+dimi-1
            for src=1:fc-1
                idx_src = find(L_symbol(tgt:n,src)) + tgt - 1;
                for rowidx=1:size(idx_src,1)
                    
                    r = idx_src(rowidx);
                    if rowidx==1 && r > tgt+dimi-1
                        break;
                    end
                    
                    L2(r,tgt) = L2(r,tgt)-L2(r,src)*transpose(L2(tgt, src));
                end
                
            end
        end
        idx = find(L_symbol(fc + dimi :n,fc)) + fc + dimi - 1;
        [blocksize , blockindex] = rowstruct2(idx, supmembership);
        
        L2(fc:fc+dimi-1,fc:fc+dimi-1) = chol(L2(fc:fc+dimi-1,fc:fc+dimi-1),'lower');
        for p = 1:size(blockindex ,1)
            r = blockindex(p);
            rb = blocksize(p);
            L2(r:r+rb-1,fc:fc+dimi-1) = L2(r:r+rb-1,fc:fc+dimi-1)/transpose(L2(fc:fc+dimi-1,fc:fc+dimi-1));
        end
    end
    
    L2 = tril(L2);
    x3=L2'\(L2\rhs);
    errorx3 = norm(x3-xtrue,'fro')/norm(xtrue,'fro')
end

for i = 1:nsuper-1
    dimi = xsuper(i+1)- xsuper(i);
    fc = xsuper(i);
    
    updater = toSupidx(find(L_symbol(fc,1:fc-1)),supmembership,xsuper);
    for r =2:dimi
        updater = union(updater,toSupidx(find(L_symbol(fc+r-1,1:fc-1)),supmembership,xsuper));
    end
    
    idx_filled = find(L_symbol(fc + dimi :n,fc)) + fc + dimi - 1;
    
    row2block = zeros(n,1);
    [blocksize , blockindex] = rowstruct2(idx_filled,supmembership);
    for b = 1:size(blockindex,1)
        frow = blockindex(b);
        rowblock = blocksize(b);
        for r = frow:frow+rowblock-1
            row2block(r) = b;
        end
    end
    
    if fullrank==0
        %compute initial ranks
        for b = 1:size(blockindex,1)
            frtarget = blockindex(b);
            rbtarget = blocksize(b);
            [u_eval, s_eval, v_eval] = svd(full(L1(frtarget:frtarget+rbtarget-1,fc:fc+dimi-1)));
            r_eval = myrank(s_eval,T);
            
            %if r_eval < dimi/3
            key = sprintf('%d_%d',i,b);
            results1(key) =  r_eval ;
            %end
        end
    end
    
    for q = 1:size(updater,1)
        k = updater(q);
        dimk = xsuper(k+1) - xsuper(k);
        j = xsuper(k);
        idx_r = find(L_symbol(fc:n,j))+fc - 1;
        [bsk, bik] = rowstruct2(idx_r,supmembership);
        
        for p = 1:size(bik,1)
            r = bik(p);
            rb = bsk(p);
            
            if fullrank==0
                btarget = row2block(r);
                [u_A, s_A, v_A] = svd(full(L1(r:r+rb-1,j:j+dimk-1)));
                r_A = myrank(s_A,T);
                k_A = u_A(:,1:r_A); l_A = s_A(1:r_A,1:r_A) ; m_A = v_A(:,1:r_A)';
                A = k_A*l_A*m_A;
                
                [u_B, s_B, v_B] = svd(full(L1(fc:fc+dimi-1, j:j+dimk-1)));
                r_B = myrank(s_B,T);
                k_B = u_B(:,1:r_B) ; l_B = s_B(1:r_B, 1:r_B) ; m_B = v_B(:,1:r_B)';
                B = k_B*l_B*m_B;
                
                D = A*transpose(B);
                [u_D, s_D, v_D] = svd(full(D));
                r_D = myrank(s_D,T);
                
                if btarget > 0
                    [u_c, s_c, v_c] = svd(full(L1(r:r+rb-1,fc:fc+dimi-1)));
                    r_c = myrank(s_c,T);
                    if r_c ~= 0
                        k_c = u_c(:,1:r_c); l_c = s_c(1:r_c,1:r_c); m_c = v_c(:,1:r_c)';
                        L1(r:r+rb-1,fc:fc+dimi -1) = (k_c *l_c * m_c) -D;
                    else
                        L1(r:r+rb-1,fc:fc+dimi -1) = L1(r:r+rb-1,fc:fc+dimi -1)-D;
                    end
                    
                    [u_C, s_C, v_C] = svd(full(L1(r:r+rb-1,fc:fc+dimi -1)));
                    r_C = myrank(s_C,T);
                    
                    
                    if ( btarget ~= 0 )
                        frtarget = blockindex(btarget);
                        rbtarget = blocksize(btarget);
                        
                        [u_eval, s_eval, v_eval] = svd(full(L1(frtarget:frtarget+rbtarget-1,fc:fc+dimi-1)));
                        r_eval = myrank(s_eval,T);
                        
                        key = sprintf('%d_%d',i,btarget);
                        results1(key) = [ results1(key); r_eval ];
                    end
                else
                    L1(r:r+rb-1,fc:fc+dimi -1) = L1(r:r+rb-1,fc:fc+dimi -1)-D;
                end
            else
                A = L1(r:r+rb-1,j:j+dimk-1);
                B = L1(fc:fc+dimi-1, j:j+dimk-1);
                D = A*transpose(B);
                L1(r:r+rb-1,fc:fc+dimi -1) = L1(r:r+rb-1,fc:fc+dimi -1)-D;
            end
        end
    end
    
    
    idx = find(L_symbol(fc + dimi :n,fc)) + fc + dimi - 1;
    [blocksize , blockindex] = rowstruct2(idx, supmembership);
    
    
    L1(fc:fc+dimi-1,fc:fc+dimi-1) = chol(L1(fc:fc+dimi-1,fc:fc+dimi-1),'lower');
    
    for p = 1:size(blockindex ,1)
        r = blockindex(p);
        rb = blocksize(p);
        if fullrank==0
            [u4, s4, v4] = svd(full(L1(r:r+rb-1,fc:fc+dimi-1)));
            r4 = myrank(s4,T);
            k4 = u4(:,1:r4); l4 = s4(1:r4,1:r4); m4 = v4(:,1:r4)';
            
            L1(r:r+rb-1,fc:fc+dimi-1) = (k4*l4*m4)/transpose(L1(fc:fc+dimi-1,fc:fc+dimi-1));
        else
            L1(r:r+rb-1,fc:fc+dimi-1) = L1(r:r+rb-1,fc:fc+dimi-1)/transpose(L1(fc:fc+dimi-1,fc:fc+dimi-1));
        end
    end
end
L1 = tril(L1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BLR implementation with rank sorting.
L = tril(Ap);
results = containers.Map('KeyType','char','ValueType','any');
for i = 1:nsuper-1
    dimi = xsuper(i+1)- xsuper(i);
    fc = xsuper(i);
    
    
    updater = toSupidx(find(L_symbol(fc,1:fc-1)),supmembership,xsuper);
    for r =2:dimi
        updater = union(updater,toSupidx(find(L_symbol(fc+r-1,1:fc-1)),supmembership,xsuper));
    end
    
    idx_filled2 = find(L_symbol(fc :n,fc)) + fc - 1;
    [blocksize2 , blockindex2] = rowstruct2(idx_filled2,supmembership);
    
    if fullrank == 0
        %compute initial ranks
        for b = 1:size(blockindex2,1)
            if b > 1
                frtarget = blockindex2(b);
                rbtarget = blocksize2(b);
                [u_eval, s_eval, v_eval] = svd(full(L(frtarget:frtarget+rbtarget-1,fc:fc+dimi-1)));
                r_eval = myrank(s_eval,T);
                
                %if r_eval < dimi/3
                key = sprintf('%d_%d',i,b-1);
                results(key) =  r_eval;
                %end
            end
        end
    end
    
    col_intervals = [];
    for s = 1:size(blockindex2,1)
        fr = blockindex2(s);
        bs = blocksize2(s);
        list = []; row_intervals = [];
        
        for t = 1:size(updater,1)
            k = updater(t);
            dimk = xsuper(k+1) - xsuper(k);
            j = xsuper(k);
            nnz = find(L_symbol(fr:fr+bs-1,j))+fr-1 ;    % point to check % you need to check column j not k which is the supernode index
            [blocksizek, blockindexk] = rowstruct2(nnz,supmembership);
            
            for m = 1:size(blockindexk,1)
                frk = blockindexk(m);
                bsk = blocksizek(m);
                lrk = frk+bsk-1;
                
                if s == 1 && lrk<fc+dimi
                    [col_intervals] = find_intervals(col_intervals,frk,lrk);
                    ncols = size(col_intervals,2);
                end
                
                [row_intervals] = find_intervals(row_intervals,frk,lrk);
                
                
                if frk >= fr && frk < fr+bs
                    for ci = 1:ncols
                        cint = col_intervals(:,ci);
                        fcc = cint(1);
                        lcc = cint(2);
                        
                        if fullrank==0
                            
                            [u_A, s_A, v_A] = svd(full(L(frk:frk+bsk-1,j:j+dimk-1)));
                            r_A = myrank(s_A,T);
                            k_A = u_A(:,1:r_A); l_A = s_A(1:r_A,1:r_A) ; m_A = v_A(:,1:r_A)';
                            A = k_A*l_A*m_A;
                            
                            
                            [u_B, s_B, v_B] = svd(full(L(fcc:lcc, j:j+dimk-1)));
                            r_B = myrank(s_B,T);
                            k_B = u_B(:,1:r_B) ; l_B = s_B(1:r_B, 1:r_B) ; m_B = v_B(:,1:r_B)';
                            B = k_B*l_B*m_B;
                            
                            D = A*transpose(B);
                            [u_D, s_D, v_D] = svd(full(D));
                            r_D = myrank(s_D,T);
                            
                            list = [list; [r_D , k ,frk, bsk,fcc,lcc]];
                        else
                            list = [list; [k , k ,frk, bsk,fcc,lcc]];
                        end
                    end
                else
                    break;
                end
                
            end
        end
        
        if size(list,1) > 0
            %divide list in multiple lists
            ncols = size(col_intervals,2);
            nrows = size(row_intervals,2);
            
            list2=cell(nrows,ncols);
            list_idx = 1;
            num_updates = size(list,1);
            
            while list_idx<=num_updates
                update = list(list_idx,:);
                ufr = update(3);
                ulr = ufr+update(4)-1;
                ufc = update(5);
                ulc = update(6);
                %find corresp row interval
                for ri = 1:nrows
                    rint = row_intervals(:,ri);
                    if ufr>=rint(1) && ulr<=rint(2)
                        break;
                    end
                end
                
                %find corresp col interval
                
                for ci = 1:ncols
                    cint = col_intervals(:,ci);
                    if ufc>=cint(1) && ulc<=cint(2)
                        break;
                    end
                end
                %list2{ri,ci} = [];
                
                update = list(list_idx,:);
                ufr = update(3);
                ulr = ufr+update(4)-1;
                ufc = update(5);
                ulc = update(6);
                
                if ufc>= cint(1) && ulc<=cint(2) && ufr>=rint(1) && ulr<=rint(2)
                    list2{ri,ci} = [ list2{ri,ci}; [update,ri,ci] ];
                end
                
                list_idx = list_idx+1;
                %list2{ri} = listrow;
            end
            update_cnt = 0;
            for ri = 1:nrows
                for ci = 1:ncols
                    listupdates = list2{ri,ci};
                    if size(listupdates,1) > 0
                        sortperm = sortrows(listupdates,[1]);
                        for m=1:size(sortperm,1)
                            tmp = sortperm(m,:);
                            r_D = sortperm(m,1);
                            k = sortperm(m,2);
                            frk = sortperm(m,3);
                            bsk = sortperm(m,4);
                            fcc = sortperm(m,5);
                            lcc = sortperm(m,6);
                            dimk = xsuper(k+1) - xsuper(k);
                            j = xsuper(k);
                            if fullrank==0
                                %recompute D
                                [u_A, s_A, v_A] = svd(full(L(frk:frk+bsk-1,j:j+dimk-1)));
                                r_A = myrank(s_A,T);
                                k_A = u_A(:,1:r_A); l_A = s_A(1:r_A,1:r_A) ; m_A = v_A(:,1:r_A)';
                                A = k_A*l_A*m_A;
                                [u_B, s_B, v_B] = svd(full(L(fcc:lcc, j:j+dimk-1)));
                                r_B = myrank(s_B,T);
                                k_B = u_B(:,1:r_B) ; l_B = s_B(1:r_B, 1:r_B) ; m_B = v_B(:,1:r_B)';
                                B = k_B*l_B*m_B;
                                D = A*transpose(B);
                                [u_D, s_D, v_D] = svd(full(D));
                                %assert(r_D == myrank(s_D,T));
                                if s > 1
                                    [u_c, s_c, v_c] = svd(full(L(frk:frk+bsk-1,fcc:lcc)));
                                    r_c = myrank(s_c,T);
                                    if r_c ~= 0
                                        k_c = u_c(:,1:r_c); l_c = s_c(1:r_c,1:r_c); m_c = v_c(:,1:r_c)';
                                        L(frk:frk+bsk-1,fcc:lcc) = (k_c *l_c * m_c) -D;
                                    else
                                        L(frk:frk+bsk-1,fcc:lcc) = L(frk:frk+bsk-1,fcc:lcc)-D;
                                    end
                                    %[u_C, s_C, v_C] = svd(full(L(frk:frk+bsk-1,fc:fc+dimi -1)));
                                    %r_C = myrank(s_C,T);
                                    
                                else
                                    L(frk:frk+bsk-1,fcc:lcc) = L(frk:frk+bsk-1,fcc:lcc)-D;
                                end
                                
                                if (  s > 1 )
                                    [u_eval, s_eval, v_eval] = svd(full(L(fr:fr+bs-1,fc:fc+dimi-1)));
                                    r_eval = myrank(s_eval,T);
                                    
                                    %if  r_eval <dimi/3
                                    key = sprintf('%d_%d',i,s-1);
                                    results(key) = [ results(key); r_eval ];
                                    %end
                                end
                            else
                                A = L(frk:frk+bsk-1,j:j+dimk-1);
                                B = L(fcc:lcc, j:j+dimk-1);
                                D = A*transpose(B);
                                L(frk:frk+bsk-1,fcc:lcc) = L(frk:frk+bsk-1,fcc:lcc)-D;
                            end
                            update_cnt = update_cnt+1;
                        end
                    end
                end
            end
            assert(update_cnt == size(list,1));
        end
    end
    
    idx = find(L_symbol(fc + dimi :n,fc)) + fc + dimi - 1;
    L(fc:fc+dimi-1,fc:fc+dimi-1) = chol(L(fc:fc+dimi-1,fc:fc+dimi-1),'lower');
    [blocksize , blockindex] = rowstruct2(idx,supmembership);
    
    for p = 1:size(blockindex ,1)
        r = blockindex(p);
        rb = blocksize(p);
        if fullrank==0
            [u4, s4, v4] = svd(full(L(r:r+rb-1,fc:fc+dimi-1)));
            r4 = myrank(s4,T);
            k4 = u4(:,1:r4); l4 = s4(1:r4,1:r4); m4 = v4(:,1:r4)';
            
            L(r:r+rb-1,fc:fc+dimi-1) = (k4*l4*m4)/transpose(L(fc:fc+dimi-1,fc:fc+dimi-1));
        else
            L(r:r+rb-1,fc:fc+dimi-1) = L(r:r+rb-1,fc:fc+dimi-1)/transpose(L(fc:fc+dimi-1,fc:fc+dimi-1));
        end
    end
end
L = tril(L);
key2 = results.keys;

x2=L_exact'\(L_exact\rhs);
x1=L1'\(L1\rhs);
x=L'\(L\rhs);

errorx2 = norm(x2-xtrue,'fro')/norm(xtrue,'fro');
errorx1 = norm(x1-xtrue,'fro')/norm(xtrue,'fro');

% errorx = norm(x-xtrue,'fro')/norm(xtrue,'fro');
% res = rhs - Ap*x;
% while norm(res,'fro') > 10e-15
%     corr=L'\(L\res);
%     x = x + corr;
%     res = rhs - Ap*x;
% end
% errorx = norm(x-xtrue,'fro')/norm(xtrue,'fro');

if fullrank==0
    totsize = 0;
    histval = [];
    for k = 1:size(key2,2)
        key = cellstr(key2(k));
        strkey = key{1};
        blkdata = sscanf(strkey,'%d_%d',2);
        supno = blkdata(1);
        blockidx=blkdata(2);
        dimsup = xsuper(supno+1)-xsuper(supno);
        %keytitle = strrep(key{1},'_',',');
        vals = cell2mat(results.values(key));
        if dimsup>1
            totsize=totsize+1;
            if(size(vals,1)>1)
                vals1 = cell2mat(results1.values(key));
                ratio = mean(vals)/mean(vals1);
                if(ratio~=1)
                    histval = [histval; ratio];
                end
            end
        end
        
    end
    
    
    h = histogram(histval)
    h.NumBins;
    [count, edges] = hist(histval);
    
    greaterthan1 = 0;
    [~, col] = find(edges >= 1);
    for t = col(1):size(edges,2)
        temp =  count(t);
        greaterthan1 = greaterthan1+temp;
    end
    
    size(histval,1)/totsize;
    
    tot2 = 0;
    for i=1:nsuper-1
        dimi = xsuper(i+1)-xsuper(i);
        if dimi>1
            tot2=tot2+1;
        end
    end
    totsize/size(key2,2);
    tot2 / (nsuper-1);
    
end

