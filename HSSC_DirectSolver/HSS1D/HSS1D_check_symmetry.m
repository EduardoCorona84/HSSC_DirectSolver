function Err = HSS1D_check_symmetry(A)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        Err = HSS1D_check_symmetry(A)
%
%    DESCRIPTION:
%        This functions finds the element-wise error of the leaves of an HSS1D structure to
%        determine if the underlying matrix is symmetric.
%
%    INPUT:
%        A is a HSS1D structure.
%
%    OUTPUT:
%        Err is the double giving the largest error of all the leaf boxes.
%                Err = max(Merr(i)) for all leaf boxes i
%            M(i) denotes the self interaction matrix of a leaf box i.
%            Merr(i) = |largest element of M'-M| / |largest element of M|
%


global BOX 

Err = zeros(size(A,2),1); 

for ibox=1:size(A,2)
    if ( (A{BOX.C1,ibox}<=0) && (A{BOX.C2,ibox}<=0) )
        Err(ibox) = max(max(abs(A{BOX.M_SELF,ibox}-A{BOX.M_SELF,ibox}.')))/max(max(abs(A{BOX.M_SELF,ibox})));
    end
end

Err = max(Err(Err>0)); 