function order = HSS1D_postorder(TREE)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        order = HSS1D_postorder(TREE)
%
%    DESCRIPTION:
%        This function outputs the box indices of HSS matrix TREE in postorder.
%
%    INPUT:
%        TREE is the {50,nboxes} cell array of a HSS matrix.
%
%    OUTPUT:
%        order is the (1 x nboxes) vector of the indices of TREE in postorder.
%


global BOX 

nboxes = size(TREE,2); 
order = zeros(1,nboxes);

ibox = 1; n = 0; 

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

%visit root
order(n+1) = 1; n = n+1;  

end

function [order,n] = visit_sub(TREE,root,order,n)
global BOX 
  
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
   
    %visit root 
   order(n+1) = root; n = n+1; 
end

