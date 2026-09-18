function [T, I] = ID(A,INPUT2)
%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%   
%   Given an input matrix A, this function determines an index 
%   vector I and a "coefficient matrix" T such that
% 
%     A(:,I) \approx A(:,I(1:k))*[I,T]
% 
%   The function has two modes of operating:
% 
%   (1) If the second input is smaller than 1, then the function
%       determines the rank k adaptively.
% 
%   (2) If the second input is larger than of equal to 1, then it
%       sets k = INPUT2
%
display(INPUT2)

if (INPUT2 < 1)
  acc = INPUT2;
  
  % QR-based rank detection. This can be switched by a more expensive but
  % accurate svd-based [T,I] = skel_col_svd(A,acc);
  [T,I] = skel_col_qr(A,acc);              
else
  k = INPUT2;
  % If k is fixed, we only need to compute the corresponding interpolation
  % operator. 
  [T,I] = skel_col_fixed(A,k);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [T, I] = skel_col_qr(A,acc)

% QR
[~, R, I] = qr(A,0);

% Using diag(R) for rank detection   
Rdg = abs(diag(R));
display(Rdg)
if ~isempty(Rdg)
Rmx = abs(Rdg(1))*acc;
k = sum(Rdg > Rmx);
T = R(1:k,1:k)\R(1:k,k+1:end);

% Power method to check interpolation error
l = min(size(A)); 
while k < l & (rand_error(A,T,I,k) > rand_norm(A,I,k)*acc)
    Rmx = 0.75*Rmx;
    k = sum(Rdg > Rmx);
    T = R(1:k,1:k)\R(1:k,k+1:end);
end

else
   T = zeros(0,size(A,2)); 
end

return

% estimate redundant matrix norm
function s = rand_norm(A,I,k)
  r = size(A,2)-k; 
  x = rand(r,1);
  x = x/norm(x);
  y = A(:,I(k+1:end))'*(A(:,I(k+1:end))*x);
  s = sqrt(norm(y));
return

% estimate approximation error
function s = rand_error(A,T,I,k)
  r = size(A,2)-k; 
  I_sk = I(1:k); I_rs = I(k+1:end); 
  
  x = rand(r,1);
  x = x/norm(x);
  y = A(:,I_sk)*(T*x) - A(:,I_rs)*x;
  y = T'*(A(:,I_sk)'*y) - A(:,I_rs)'*y;
  s = sqrt(norm(y));
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [T, I] = skel_col_svd(A,acc)

%QR
[~, R, I] = qr(A,0);

% 'Expensive / accurate' rank detection way via SVD of R
ss = svd(A);
k  = sum(ss > acc);
[U,D,V] = svd(R(1:k,1:k));
q = sum(diag(D) > acc);
T = V(:,1:q)*(D(1:q,1:q)\(U(:,1:q)'))*R(1:k, (k+1):size(R,2));

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [T, I] = skel_col_fixed(A,k)

[~, R, I] = qr(A,0);
[U,D,V] = svd(R(1:k,1:k));
q = sum(diag(D) > 1e-12);
T = V(:,1:q)*(D(1:q,1:q)\(U(:,1:q)'))*R(1:k, (k+1):size(R,2));

return
