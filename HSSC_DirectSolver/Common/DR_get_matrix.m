function [A,A_dr] = DR_get_matrix(xx,h,n,rk,order)
%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%

  
ntot = size(xx,2);

%%% Set up the "regular weights".
%wwvec = [0.5*h,h*ones(1,n-2),0.5*h];
wwvec = h*ones(1,n); 
WW    = wwvec' * wwvec;
ww    = WW(:)';

%%% Compute the DR corrections.
[J,dr_weights,iistencil] = DR_stencil(h,rk,n,order);

%%% Compute the "regular" part of the matrix.
DD1            = xx(1,:)'*ones(1,ntot) - ones(ntot,1)*xx(1,:);
DD2            = xx(2,:)'*ones(1,ntot) - ones(ntot,1)*xx(2,:);
DD             = sqrt(DD1.*DD1 + DD2.*DD2);

size(DD)
size(ww)

A_reg          = besselh(0,rk*DD + eye(ntot)).*(ones(ntot,1)*ww);
inddiag        = (1:ntot) + ntot*(0:(ntot-1));
A_reg(inddiag) = zeros(1,ntot);

%%% Compute the "correction" matrix.
fprintf(1,'Computing the "correction weight matrix".\n')
border = max(max(abs(iistencil(1:2,:))));
INDINT = n*ones(n-2*border,1)*(border:(n-border-1)) + ...
         ((border+1):(n-border))' * ones(1,n-2*border);
indint = INDINT(:)';
II1    = ones(length(J),1)*indint;
II2    = J'               *ones(1,length(indint)) + II1;
DDD    = (dr_weights.')   *ones(1,length(indint));
A_dr   = sparse(II1(:),II2(:),DDD(:),ntot,ntot);

%%% Assemble the matrix
A = A_reg + full(A_dr);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

