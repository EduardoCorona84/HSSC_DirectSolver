function [HSSinv,FTinv] = HSS1D_invertp_nsym(K_HSS,acc)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        [HSSinv,FTinv] = HSS1D_invertp_nsym(K_HSS,acc)
%
%    DESCRIPTION:
%        This function inverts a non-symmetric HSS1D matrix K_HSS to desired accuracy acc using
%        pseudo inverse to invert blocks in the structure.
%
%    INPUT:
%        K_HSS is the {50 x nboxes} cell array that stores the HSS matrix to be inverted.
%        acc (double) is the desired accuracy
%
%    OUTPUT:
%        HSSinv is a {3 x nboxes} cell array, corresponding to the HSSinv structure in 
%        the Martinsson-Rokhlin paper.  
% 
%        HSSinv{BOX.Linv,BOX.Rinv and BOX.Dinv,nboxes) correspond to the diagonal
%        blocks of the telescopic factorization of the inverse as explained in the
%        CMZ_2012 paper. That is, at each level the inverse can be written as: 
%
%                    (K_HSS)^{-1} = Dinv + Linv*(M+E)^{-1}*Rinv
%
%        FTinv is the schur complement matrix of the HSSinv structure at the 
%        top resulting from the merge of its children's interactions. 
%        It can also be considered Dinv for the top box. 
%

% Cell index variable names
global BOX 
% Extra names for inverse cell HSSinv
BOX.Linv = 1; 
BOX.Rinv = 2; 
BOX.Dinv = 3; 

nboxes  = size(K_HSS,2);
HSSinv     = cell(3,nboxes);
E_VEC = cell(1,nboxes);

% Loop over all K_HSS, from finest to coarser.
for ibox = nboxes:(-1):2
    
  % Assemble the diagonal matrix.
  if ((K_HSS{BOX.C1,ibox} <= 0) && (K_HSS{BOX.C2,ibox} <= 0)) % ibox is a leaf.
    F    = K_HSS{BOX.M_SELF,ibox};
  elseif ((K_HSS{BOX.C1,ibox}>0) && (K_HSS{BOX.C2,ibox}<=0)) % ibox has one son - the left one.
    ison1 = K_HSS{BOX.C1,ibox};
    F    = E_VEC{ison1};
  elseif ((K_HSS{BOX.C1,ibox}<=0) && (K_HSS{BOX.C2,ibox}>0)) % ibox has one son - the right one.
    ison2 = K_HSS{BOX.C2,ibox};
    F    = E_VEC{ison2};
  elseif ((K_HSS{BOX.C1,ibox}>0) && (K_HSS{BOX.C2,ibox}>0)) % ibox has two sons.
    ison1 = K_HSS{BOX.C1,ibox};
    ison2 = K_HSS{BOX.C2,ibox};
    F    = [E_VEC{ison1},K_HSS{BOX.M_SIB,ison1};...
             K_HSS{BOX.M_SIB,ison2},E_VEC{ison2}];
  end
  
  % Construct the interpolation matrices L and R from K_HSS. 
  n                  = K_HSS{BOX.NSKEL,ibox};
  k                  = K_HSS{BOX.KSKEL,ibox};
  Jout               = K_HSS{BOX.J_UP,ibox};
  Jin                = K_HSS{BOX.J_DN,ibox};
  R                  = zeros(k,n);
  R(:,    Jout(1:k)) = eye(k);
  R(:,Jout((k+1):n)) = K_HSS{BOX.T_UP,ibox};
  L                  = zeros(n,k);
  L(Jin(1:k),    :)  = eye(k);
  L(Jin((k+1):n),:)  = K_HSS{BOX.T_DN,ibox}.';
  
  % Construct the various projection maps.
  Finv            = pinv(F,acc);
  EI             = R*(F\L);
  E              = pinv(EI,acc);
  E_VEC{ibox}    = E;
  HSSinv{BOX.Linv,ibox}      = F\(L*E);
  HSSinv{BOX.Rinv,ibox}      = (EI\R)*Finv;
  HSSinv{BOX.Dinv,ibox}      = Finv - HSSinv{BOX.Linv,ibox}*(R*Finv);

end

% Assemble the "top matrix" and invert it:
FT = [E_VEC{2},K_HSS{BOX.M_SIB,2};...
      K_HSS{BOX.M_SIB,3},E_VEC{3}];

FTinv = pinv(FT,acc);