function A_NEW = HSS1D_symmetrize_fsym(A)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        A_NEW = HSS1D_symmetrize_fsym(A)
%
%    DESCRIPTION:
%        This function enforces that sibling and self interactions are symmetric.
%
%    INPUT:
%        A is the {50 x nbox} cell array that encodes the HSS info of a matrix.
%
%    OUTPUT:
%        A_NEW is the {50 x nbox} cell array that encodes the symmetrized matrix.
%


global BOX 

A_NEW = A; 

for ibox=1:size(A,2)
    if ( (A{BOX.C1,ibox}>0) && (A{BOX.C2,ibox}>0) )
    ison1 = A{BOX.C1,ibox};
    ison2 = A{BOX.C2,ibox};
    A_NEW{BOX.M_SIB,ison2} = A{BOX.M_SIB,ison1}.';
    end
    if ( (A{BOX.C1,ibox}<=0) && (A{BOX.C2,ibox}<=0) )
        A_NEW{BOX.M_SELF,ibox} = (A{BOX.M_SELF,ibox}+A{BOX.M_SELF,ibox}.')/2;
    end
end
