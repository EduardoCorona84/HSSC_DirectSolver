function lHSS = levHSS(TREE,n_cut,trinv)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%     FUNCTION CALL:
%         lhss = levHSS(TREE,n_cut,trinv)
%
%     DESCRIPTION:
%         This function does a downward pass of the tree and finds the first level that has a
%         box with less than n_cut points. If no box has more than n_cut points, then lHSS = -1.
%         The return value determines the cutoff level used in the O(N) version of the COMPRESS
%         and INVERT stages; matrices in boxes above level lhss are compressed, and those below
%         or on level lhss are stored densely. This function needs to be called by the user for
%         the O(N) algorithm and the return value needs to be stored in params.lev_hss
%         (the struct created in HSS_tree_parameters.m).
%
%     INPUT:
%         TREE    <struct>    HSS binary tree data (e.g. output of HSS2D_bintree_*).
%         n_cut   <int>       Threshold box size (e.g. params.n_cut from HSS_tree_parameters).
%         triinv  <bool>      Tree traversal per box or per level (for translation invariance).
%
%     OUTPUT:
%         lhss    <int>       First level in the tree that doesn't use HSS1D speedup.
%





lev=0;
if trinv == 0
    n_min = max([TREE.BOX(1).k,TREE.BOX(2).k,TREE.BOX(3).k]);   
else
    n_min = TREE.LEV(1).k; 
end

while (n_min>n_cut && lev<TREE.depth - 1)
    lev = lev+1;
    k_sk = zeros(TREE.numlev(lev+1),1);
    
    if trinv==0
        for i=1:TREE.numlev(lev+1)
            Nbox = TREE.box_numbers(lev+1,i);
            k_sk(i) = TREE.BOX(Nbox).k;
        end
    
        n_min = min(k_sk); 
    else
        n_min = TREE.LEV(lev+1).k;
    end
        
end

if lev>0
    if mod(lev,2)==1
        lHSS = lev-1;
    else
        lHSS = lev;
    end
else
    lHSS = -1; 
end
