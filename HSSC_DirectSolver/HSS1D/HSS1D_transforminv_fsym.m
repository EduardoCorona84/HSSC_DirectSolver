function Kinv_HSS = HSS1D_transforminv_fsym(K_HSS,HSSinv,FTinv)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2008-09 P.G. Martinsson (C) 2011-2013 E. Corona, P.G. Martinsson, D. Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        Kinv_HSS = HSS1D_transforminv_fsym(K_HSS,HSSinv,FTinv)
%
%    DESCRIPTION:
%        This function converts the inverse of a full-symmetric HSS matrix A in the HSSinv format to
%        the same HSS format as A. This is performed by computing sibling interactions for non-leaf
%        boxes and converting Linv and Rinv matrices into interpolation matrices that look like [I T].
%
%    INPUT:
%        K_HSS is the {50 x nboxes} cell array that encodes the HSS structure of A.
%        HSSinv, FTinv encode the inverse of A in the HSSinv format and are the output of function
%            HSS1D_invert_fsym(..).
%
%    OUTPUT:
%        Kinv_HSS is the {50 x nboxes} cell array that encodes the HSS structure of inv(A).
%

% Cell index variable names
global BOX 
% Extra names for inverse cell HSSinv
BOX.Linv = 1; 
BOX.Rinv = 2; 
BOX.Dinv = 3;

nboxes = size(K_HSS,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Move all "diagonal" contributions down to Dinv at the lowest level. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Move the diagonal from the smallest "full" matrix to Dinv at level 1:
ison1 = 2;
ison2 = 3;
k1    = K_HSS{BOX.KSKEL,ison1};
k2    = K_HSS{BOX.KSKEL,ison2};
ind1  =      (1:k1);
ind2  = k1 + (1:k2);
HSSinv{BOX.Dinv,ison1} = HSSinv{BOX.Dinv,ison1} + HSSinv{BOX.Rinv,ison1}.'*FTinv(ind1,ind1)*HSSinv{BOX.Rinv,ison1};
HSSinv{BOX.Dinv,ison2} = HSSinv{BOX.Dinv,ison2} + HSSinv{BOX.Rinv,ison2}.'*FTinv(ind2,ind2)*HSSinv{BOX.Rinv,ison2};
FTinv(ind1,ind1) = zeros(k1,k1);
FTinv(ind2,ind2) = zeros(k2,k2);
% Continue moving the diagonal blocks down towards finer levels.
for ibox = 2:nboxes
  if ( (K_HSS{BOX.C1,ibox} >  0) & (K_HSS{BOX.C2,ibox} >  0)) % ibox has two sons
    ison1 = K_HSS{BOX.C1,ibox};
    ison2 = K_HSS{BOX.C2,ibox};
    k1    = K_HSS{BOX.KSKEL,ison1};
    k2    = K_HSS{BOX.KSKEL,ison2};
    ind1  =      (1:k1);
    ind2  = k1 + (1:k2);
    HSSinv{BOX.Dinv,ison1} = HSSinv{BOX.Dinv,ison1} + HSSinv{BOX.Rinv,ison1}.'*HSSinv{BOX.Dinv,ibox}(ind1,ind1)*HSSinv{BOX.Rinv,ison1};
    HSSinv{BOX.Dinv,ison2} = HSSinv{BOX.Dinv,ison2} + HSSinv{BOX.Rinv,ison2}.'*HSSinv{BOX.Dinv,ibox}(ind2,ind2)*HSSinv{BOX.Rinv,ison2};
    HSSinv{BOX.Dinv,ibox}(ind1,ind1) = zeros(k1,k1);
    HSSinv{BOX.Dinv,ibox}(ind2,ind2) = zeros(k2,k2);
  elseif ( (K_HSS{BOX.C1,ibox} >  0) & (K_HSS{BOX.C2,ibox} <= 0)) % ibox has one son - the left.
    ison1        = K_HSS{BOX.C1,ibox};
    k1           = K_HSS{BOX.KSKEL,ison1};
    HSSinv{BOX.Dinv,ison1} = HSSinv{BOX.Dinv,ison1} + HSSinv{BOX.Rinv,ison1}.'*HSSinv{BOX.Dinv,ibox}*HSSinv{BOX.Rinv,ison1};
    HSSinv{BOX.Dinv,ibox}  = zeros(k1,k1);
  elseif ( (K_HSS{BOX.C1,ibox} <= 0) & (K_HSS{BOX.C2,ibox} >  0)) % ibox has one son - the right.
    ison2        = K_HSS{BOX.C2,ibox};
    k2           = K_HSS{BOX.KSKEL,ison2};
    HSSinv{BOX.Dinv,ison2} = HSSinv{BOX.Dinv,ison2} + HSSinv{BOX.Rinv,ison2}.'*HSSinv{BOX.Dinv,ibox}*HSSinv{BOX.Rinv,ison2};
    HSSinv{BOX.Dinv,ibox}  = zeros(k2,k2);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Transform the factorization to make the "E" and "F" matrices %%%
%%% "skeletonization matrices".                                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ibox = nboxes:(-1):2
  n  = K_HSS{BOX.NSKEL,ibox};
  k  = K_HSS{BOX.KSKEL,ibox};
  J  = K_HSS{BOX.J_UP,ibox};
  R  = HSSinv{BOX.Rinv,ibox};
  Q  = R(:,J(1:k));
  Rnew = zeros(k,n);
  Rnew(:,J(1:k    )) = eye(k);
  Rnew(:,J((k+1):n)) = Q\R(:,J((k+1):n));
  HSSinv{BOX.Rinv,ibox} = Rnew;
  ifath = K_HSS{BOX.PARENT,ibox};
  % Check whether ibox is the left or the right son of its father.
  if ( (ibox==K_HSS{BOX.C2,ifath}) & (K_HSS{BOX.C1,ifath}>0)) % ibox is the right son AND there is a left son.
    ison1fath = K_HSS{BOX.C1,ifath};
    ind       = K_HSS{BOX.KSKEL,ison1fath} + (1:k);
  else
    ind       = 1:k;
  end
  % Check whether we are at level 1 yet.
  % If not, then update E,F,G on the next coarser level.
  % If yes, then update FTinv.
  if (K_HSS{BOX.LEVEL,ibox} > 1)
    HSSinv{BOX.Rinv,ifath}(:,ind) =    HSSinv{BOX.Rinv,ifath}(:,ind)*Q;
    HSSinv{BOX.Dinv,ifath}(ind,:) = Q.'*HSSinv{BOX.Dinv,ifath}(ind,:)  ;
    HSSinv{BOX.Dinv,ifath}(:,ind) =    HSSinv{BOX.Dinv,ifath}(:,ind)*Q;
  else % ibox is at level 1.
    FTinv(ind,:) = Q.'*FTinv(ind,:)  ;
    FTinv(:,ind) =    FTinv(:,ind)*Q;
  end
end

% Transfer the information from E,F,G,FTinv to a new Hudson structure.
Kinv_HSS = K_HSS;

% Extract the "T" matrices.
for ibox = 2:nboxes
  n = K_HSS{BOX.NSKEL,ibox};
  k = K_HSS{BOX.KSKEL,ibox};
  J = K_HSS{BOX.J_UP,ibox};
  R = HSSinv{BOX.Rinv,ibox};
  T = R(:,J((k+1):n));
  Kinv_HSS{BOX.T_UP,ibox} = T;
  Kinv_HSS{BOX.J_UP,ibox} = J;
end

% Extract the sibling interaction matrices on level 1 from FTinv.
ison1 = 2;
ison2 = 3;
ind1  = 1:K_HSS{BOX.KSKEL,ison1};
ind2  = K_HSS{BOX.KSKEL,ison1} + (1:K_HSS{BOX.KSKEL,ison2});
Kinv_HSS{BOX.M_SIB,ison1} = FTinv(ind1,ind2);
Kinv_HSS{BOX.M_SIB,ison2} = Kinv_HSS{BOX.M_SIB,ison1}.';
% Extract information from the "G" matrices for all boxes.
% For non-leaf K_HSS, this information becomes sibling interaction matrices.
% For leah K_HSS, the G matrix becomes the matrix of self interaction.
for ibox = 2:nboxes
  if ( (K_HSS{BOX.C1,ibox} > 0) & (K_HSS{BOX.C2,ibox} > 0) )
    ison1 = K_HSS{BOX.C1,ibox};
    ison2 = K_HSS{BOX.C2,ibox};
    ind1  = 1:K_HSS{BOX.KSKEL,ison1};
    ind2  = K_HSS{BOX.KSKEL,ison1} + (1:K_HSS{BOX.KSKEL,ison2});
    Kinv_HSS{BOX.M_SIB,ison1} = HSSinv{BOX.Dinv,ibox}(ind1,ind2);
    Kinv_HSS{BOX.M_SIB,ison2} = Kinv_HSS{BOX.M_SIB,ison1}.';
  elseif ( (K_HSS{BOX.C1,ibox}<=0) & (K_HSS{BOX.C2,ibox}<=0) )
    Kinv_HSS{BOX.M_SELF,ibox} = HSSinv{BOX.Dinv,ibox};
    %Kinv_HSS{BOX.M_SELF,ibox}  = (HSSinv{BOX.Dinv,ibox}+HSSinv{BOX.Dinv,ibox}')/2;
  end
end