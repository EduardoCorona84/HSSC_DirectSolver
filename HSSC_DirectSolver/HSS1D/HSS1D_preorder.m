function order = HSS1D_preorder(TREE)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        order = HSS1D_preorder(TREE)
%
%    DESCRIPTION:
%        This function outputs the box indices of the HSS matrix TREE in preorder.
%
%    INPUT:
%        TREE is the {50 x nboxes} cell array of a HSS matrix.
%
%    OUTPUT:
%        order is the (1 x nboxes) vector of the indices of TREE in preorder.
%


global BOX 

nboxes = size(TREE,2); 
order = zeros(1,nboxes);

%visit root
order(1) = 1; n = 1; ibox = 1; 

%visit left subtree
ison1 = TREE{BOX.C1,ibox};
if (ison1 >0)
[order,n] = visit_sub(TREE,ison1,order,n); 
end

%visit right subtree
ison2 = TREE{BOX.C2,ibox}; 
if (ison2 >0)
[order,n] = visit_sub(TREE,ison2,order,n); 
end

end

function [order,n] = visit_sub(TREE,root,order,n)
global BOX 
   %visit root 
   order(n+1) = root; 
   n = n+1; 
   ison1 = TREE{BOX.C1,root};
   
   %visit left subtree
   if (ison1>0)
       [order,n] = visit_sub(TREE,ison1,order,n); 
   end
   
   %visit right subtree
   ison2 = TREE{BOX.C2,root};
   if (ison2>0)
       [order,n] = visit_sub(TREE,ison2,order,n); 
   end
end

