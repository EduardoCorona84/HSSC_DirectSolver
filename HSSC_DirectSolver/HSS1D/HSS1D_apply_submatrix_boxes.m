function [Box_dn,Box_up] = HSS1D_apply_submatrix_boxes(TREE,I_dn,I_up)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        [Box_dn,Box_up] = HSS1D_apply_submatrix_boxes(TREE,I_dn,I_up)
%
%    DESCRIPTION:
%        This function finds the subtrees in the HSS1D structure TREE that correspond to the
%        submatrix A[I_dn,I_up].
%
%    INPUT:
%        A denotes the matrix encoded in HSS1D structure TREE.
%        I_dn / I_up are boolean vectors that specify the row / column index sets of the submatrix
%            of A.
%
%    OUTPUT:
%        Box_dn / Box_up are the outputted boolean index vectors that specify the boxes of TREE
%            corresponding to submatrix A[I_dn,I_up].
%


global BOX 

nboxes = size(TREE,2); 

Box_dn = false(nboxes,1); 
Box_up = Box_dn; 

for ibox = nboxes:-1:1
    if ( (TREE{BOX.C1,ibox}<=0) && (TREE{BOX.C2,ibox}<=0) )
       Idnbox = I_dn((TREE{BOX.END1,ibox} - 1 + (1:TREE{BOX.END2,ibox})));
       Iupbox = I_up((TREE{BOX.END1,ibox} - 1 + (1:TREE{BOX.END2,ibox}))); 
       
       N1dn = sum(Idnbox); N1up = sum(Iupbox); 
       if N1dn>0
           Box_dn(ibox) = true; 
       end
       if N1up>0
           Box_up(ibox) = true; 
       end
    else
        ison1 = TREE{BOX.C1,ibox}; 
        ison2 = TREE{BOX.C2,ibox}; 
        if ison1>0 && ison2>0
            if (Box_dn(ison1)==true || Box_dn(ison2)==true)
                Box_dn(ibox)=true;
            end
            if (Box_up(ison1)==true || Box_up(ison2)==true)
                Box_up(ibox)=true;
            end
        else
           ison = max(ison1,ison2); 
           if (Box_dn(ison)==true)
              Box_dn(ibox)=true;  
           end
           if (Box_up(ison)==true)
              Box_up(ibox)=true;  
           end
        end
    end
end