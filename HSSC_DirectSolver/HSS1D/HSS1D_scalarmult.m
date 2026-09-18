function A = HSS1D_scalarmult(c,A)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        A = HSS1D_scalarmult(c,A)
%
%    DESCRIPTION:
%        This function performs a matrix-scalar product where the input and matrix are in HSS1D
%        form : A = cA
%
%    INPUT:
%        c (double) scalar
%        A is the cell array containing an HSS matrix.
%
%    OUTPUT:
%        A is the cell array containing HSS matrix cA.
%


global BOX 

for ibox=1:size(A,2)
    if ( (A{BOX.C1,ibox}>0) && (A{BOX.C2,ibox}>0) )
    ison1 = A{BOX.C1,ibox};
    ison2 = A{BOX.C2,ibox};
    A{BOX.M_SIB,ison1} = c*A{BOX.M_SIB,ison1};
    A{BOX.M_SIB,ison2} = c*A{BOX.M_SIB,ison2};
    end
    if ( (A{BOX.C1,ibox}<=0) && (A{BOX.C2,ibox}<=0) )
        A{BOX.M_SELF,ibox} = c*A{BOX.M_SELF,ibox};
    end
end