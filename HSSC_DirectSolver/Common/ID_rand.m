function [T, I] = ID_rand(m,n,matvecA,matvecAt,dist,params,k0,C)
%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
% 		This algorithm is a variation of the standard ID code which uses randomized
% 		methods to compute an interpolative decomposition, that is: 
%		
% 		  A(:,I) \approx A(:,I(1:k))*[I,T]
%		
% 		The function has two modes of operating:
%		
% 		(1) If the second input is smaller than 1, then the function
% 		    determines the rank k adaptively via randomized methods.
%		
% 		(2) If the second input is larger than of equal to 1, then it
% 		    sets k = params
%		
% 		The structure of this code is based on the iddecomp code (Martinsson,
% 		OMNI). It uses the RandQZ and double GS functions. 
% 		Programmer: Eduardo Corona
% 		Numerical Linear Algebra with Probability, Mark Tygert 2010
%

if (nargin == 8)
  tol = params;
  [Z,k] = Rand_AT_sample(m,n,matvecA,matvecAt,dist,tol,k0,C);
  [T,I] = skel_colk(Z,k);
else
  k = params;
  [Z,~] = Rand_AT_sample(m,n,matvecA,matvecAt,dist,k);
  k = size(Z,1); 
  [T,I] = skel_colk(Z,k);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [T, I] = skel_colk(A,k)

[~, R, I] = qr(A,0);

[U,D,V] = svd(R(1:k,1:k));
q = sum(diag(D) > 1e-12);
T = V(:,1:q)*(D(1:q,1:q)\(U(:,1:q)'))*R(1:k, (k+1):size(R,2));

return
