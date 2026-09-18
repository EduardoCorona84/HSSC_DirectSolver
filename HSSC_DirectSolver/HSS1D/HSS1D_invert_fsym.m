function [HSSinv,FTinv] = HSS1D_invert_fsym(K_HSS,acc)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2008-09 P.G. Martinsson (C) 2011-2013 E. Corona, P.G. Martinsson, D. Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        [HSSinv,FTinv] = HSS1D_invert_fsym(K_HSS,acc)
%
%    DESCRIPTION:
%        This function inverts a fully-symmetric HSS1D matrix K_HSS to desired accuracy acc.
%
%    INPUT:
%        K_HSS is the {50 x nboxes} cell array that stores the HSS matrix to be inverted.
%        acc (double) is the desired accuracy.
%
%    OUTPUT:
%        HSSinv is a {3 x nboxes} cell array, corresponding to the EFG structure in 
%        the Martinsson-Rokhlin paper.  
% 
%        HSSinv{BOX.Linv,BOX.Rinv and BOX.Dinv,nboxes) correspond to the diagonal
%        blocks of the telescopic factorization of the inverse as explained in the
%        CMZ_2012 paper. Linv matrices are not explicitly computed because
%            of symmetry (Linv = Rinv.'). That is, at each level the
%            inverse can be written as: 
%
%                    (K_HSS)^{-1} = Dinv + Linv*(M+E)^{-1}*Rinv
%
%        FTinv is the schur complement matrix of the HSSinv structure at the top resulting from the merge of its
%            children's interactions. It can also be considered Dinv for
%            the top box. 
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
  if ((K_HSS{BOX.C1,ibox} <= 0) & (K_HSS{BOX.C2,ibox} <= 0)) % ibox is a leaf.
    F    = K_HSS{BOX.M_SELF,ibox};
  elseif ((K_HSS{BOX.C1,ibox}>0) & (K_HSS{BOX.C2,ibox}<=0)) % ibox has one son - the left one.
    ison1 = K_HSS{BOX.C1,ibox};
    F    = E_VEC{ison1};
  elseif ((K_HSS{BOX.C1,ibox}<=0) & (K_HSS{BOX.C2,ibox}>0)) % ibox has one son - the right one.
    ison2 = K_HSS{BOX.C2,ibox};
    F    = E_VEC{ison2};
  elseif ((K_HSS{BOX.C1,ibox}>0) & (K_HSS{BOX.C2,ibox}>0)) % ibox has two sons.
    ison1 = K_HSS{BOX.C1,ibox};
    ison2 = K_HSS{BOX.C2,ibox};
    F    = [E_VEC{ison1},K_HSS{BOX.M_SIB,ison1};...
             K_HSS{BOX.M_SIB,ison1}.',E_VEC{ison2}];
  end
  
  % Construct the interpolation L matrices from K_HSS. Since the kernel is symmetric, R
  % is the same as L.'. 
  n                  = K_HSS{BOX.NSKEL,ibox};
  k                  = K_HSS{BOX.KSKEL,ibox};
  Jout               = K_HSS{BOX.J_UP,ibox};
  L                  = zeros(n,k);
  L(Jout(1:k),    :) = eye(k);
  L(Jout((k+1):n),:) = K_HSS{BOX.T_UP,ibox}.';
  
  % Construct scattering / schur complement matrices E and Finv, and use them
  % to build the factors Linv, Rinv and Dinv for the inverse telescopic
  % factorization. 
  Finv            = inv(F);
  E              = inv(L.'*Finv*L);
  E_VEC{ibox}    = Einv;
  % HSSinv factors. Again because of symmetry, Linv is Rinv.'. 
  HSSinv{BOX.Rinv,ibox}      = E*(L.')*Finv;
  HSSinv{BOX.Dinv,ibox}      = Finv - (Finv*L)*HSSinv{BOX.Rinv,ibox};

end

% Assemble the "top F matrix" and invert it:
FT = [E_VEC{2},  K_HSS{BOX.M_SIB,2};...
      K_HSS{BOX.M_SIB,2}.',E_VEC{3}];

FTinv = inv(FT);
