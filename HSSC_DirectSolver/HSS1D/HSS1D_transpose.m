function At_HSS = HSS1D_transpose(A_HSS)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        At_HSS = HSS1D_transpose(A_HSS)
%
%    DESCRIPTION:
%        This function returns the HSS structure of the transpose of matrix A.
%
%    INPUT:
%        A_HSS is the {50 x nbox} cell array that encodes the HSS structure of matrix A.
%
%    OUTPUT:
%        At_HSS is the {50 x nbox} cell array that encodes the HSS structure of matrix A.'.
%


global BOX 

% Initialize with A_HSS of A
At_HSS = A_HSS;
nboxes = size(A_HSS,2);

% Ingoing and Outgoing trees are the same

% Exchange ingoing and outgoing skeleton information

At_HSS(BOX.I_SKUP,:)=A_HSS(BOX.I_SKDN,:);
At_HSS(21,:)=A_HSS(23,:);
At_HSS(BOX.I_SKDN,:)=A_HSS(BOX.I_SKUP,:);
At_HSS(23,:)=A_HSS(21,:);


% Exchange and transpose matrices

At_HSS(BOX.T_UP,:) = A_HSS(BOX.T_DN,:);
At_HSS(BOX.J_UP,:) = A_HSS(BOX.J_DN,:);
At_HSS(BOX.T_DN,:) = A_HSS(BOX.T_UP,:);
At_HSS(BOX.J_DN,:) = A_HSS(BOX.J_UP,:);

for ibox=1:nboxes
    if ( (A_HSS{BOX.C1,ibox}>0) & (A_HSS{BOX.C2,ibox}>0) )
    ison1 = A_HSS{BOX.C1,ibox};
    ison2 = A_HSS{BOX.C2,ibox};
    At_HSS{BOX.M_SIB,ison1} = A_HSS{BOX.M_SIB,ison2}.';
    At_HSS{BOX.M_SIB,ison2} = A_HSS{BOX.M_SIB,ison1}.';
    end
    At_HSS{BOX.M_SELF,ibox} = A_HSS{BOX.M_SELF,ibox}.';
end
