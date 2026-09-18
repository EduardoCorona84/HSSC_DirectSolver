function xx = HSS1D_applyinv_nsym(K_HSS,HSSinv,FTinv,ff,dim)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2008-09 P.G. Martinsson (C) 2011-2013 E. Corona, P.G. Martinsson, D. Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        xx = HSS1D_applyinv_nsym(K_HSS,HSSinv,FTinv,ff,dim)
%
%    DESCRIPTION:
%        This function uses the HSSinv structure {HSSinv,FTinv} of a non-symmetric HSS1D matrix A to
%        apply its inverse : xx = inv(A)ff.
%
%    INPUT:
%        A is the matrix encoded by HSS1D structure K_HSS
%        HSSinv / FTinv are the hierarchical components of inv(A). They can be precomputed by
%            [HSSinv, FTinv] = HSS1D_invert_nsym(K_HSS,acc)
%        ff is the vector to be left-multiplied by inv(A)
%        dim (int) is the dimension of the kernal.
%
%    OUTPUT:
%        xx is the vector resulting from the matrix-vector multiply
%            xx = inv(A)ff
%

% Cell index variable names
global BOX 
% Extra names for inverse cell HSSinv
BOX.Linv = 1; 
BOX.Rinv = 2; 
BOX.Dinv = 3; 

if nargin == 4
    dim = 1; 
end

nboxes = size(K_HSS,2);
FIELDS = cell(2,nboxes);

% The upwards pass - collect all "outgoing fields". 
% This corresponds to multiplication by Rinv matrices. 
for ibox = nboxes:(-1):2
    
  if ( (K_HSS{BOX.C1,ibox} <= 0) && (K_HSS{BOX.C2,ibox} <= 0)) % ibox is a leaf
    ind = dim*K_HSS{BOX.END1,ibox} - dim + (1:dim*K_HSS{BOX.END2,ibox});
    FIELDS{1,ibox} = HSSinv{BOX.Rinv,ibox}*ff(ind,:);
    
  elseif ( (K_HSS{BOX.C1,ibox} > 0) && (K_HSS{BOX.C2,ibox} > 0))  % ibox has two sons.
    ison1       = K_HSS{BOX.C1,ibox};
    ison2       = K_HSS{BOX.C2,ibox};
    field_in    = [FIELDS{1,ison1};FIELDS{1,ison2}];
    FIELDS{1,ibox} = HSSinv{BOX.Rinv,ibox}*field_in;
  else                                                  % ibox has one son
      ison = max(K_HSS{BOX.C1,ibox},K_HSS{BOX.C2,ibox});
      FIELDS{1,ibox} = HSSinv{BOX.Rinv,ibox}*FIELDS{1,ison};
  end
end

% Process the top-level.
f_top = [FIELDS{1,2};FIELDS{1,3}];
x_top = FTinv*f_top;
FIELDS{2,2} = x_top(1:size(HSSinv{BOX.Linv,2},2),:);
FIELDS{2,3} = x_top((size(HSSinv{BOX.Linv,2},2)+1):end,:);

% Downwards pass: collection of "ingoing fields" in FIELDS{2,ibox}. 
% This involves operations of the form: Linv*FIELDS{2,ibox}+Dinv*[FIELDS{1,isons}]
xx = zeros(size(ff));

for ibox = 2:nboxes
  if ( (K_HSS{BOX.C1,ibox} <= 0) && (K_HSS{BOX.C2,ibox} <= 0)) % ibox is a leaf
    ind = dim*K_HSS{BOX.END1,ibox} - dim + (1:dim*K_HSS{BOX.END2,ibox});
    xx(ind,:) = HSSinv{BOX.Linv,ibox}*FIELDS{2,ibox} + HSSinv{BOX.Dinv,ibox}*ff(ind,:);
  elseif ( (K_HSS{BOX.C1,ibox} > 0) && (K_HSS{BOX.C2,ibox} > 0))  % ibox has two sons.
    ison1  = K_HSS{BOX.C1,ibox};
    ison2  = K_HSS{BOX.C2,ibox};
    k1     = K_HSS{BOX.KSKEL,ison1};
    k2     = K_HSS{BOX.KSKEL,ison2};
    x_long = HSSinv{BOX.Linv,ibox}*FIELDS{2,ibox} + HSSinv{BOX.Dinv,ibox}*[FIELDS{1,ison1};FIELDS{1,ison2}];
    FIELDS{2,ison1} = x_long(     (1:k1),:);
    FIELDS{2,ison2} = x_long(k1 + (1:k2),:);
  else
      ison = max(K_HSS{BOX.C1,ibox},K_HSS{BOX.C2,ibox});
      FIELDS{2,ison} = HSSinv{BOX.Linv,ibox}*FIELDS{2,ibox} + HSSinv{BOX.Dinv,ibox}*FIELDS{1,ison};
  end
end
