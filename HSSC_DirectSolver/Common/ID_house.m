function [rss,vn,scal] = ID_house(n,x)
%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%     
%     Constructs the vector vn with vn(1) = 1, and the scalar scal, 
%     such that the obviously self-adjoint
%     H := identity_matrix - scal * vn * adjoint(vn) is unitary,
%     the absolute value of the first entry of Hx
%     is the root-sum-square of the entries of x,
%     and all other entries of Hx are zero
%     (H is the Householder matrix corresponding to x).
%     
%     Input:
%       x -- vector to reflect into its first component
%     
%     Output:
%       css -- root-sum-square of the entries of x * the phase of x(1)
%       vn -- entries 2 to n of the Householder vector vn;
%       vn(1) is assumed to be 1
%       scal -- scalar multiplying vn * adjoint(vn);
%       
%       scal = 2/(1 + |vn(2)|^2 + ... + |vn(n)|^2) when vn(2), ..., vn(n) don't all vanish;
%       
%       scal = 0 when vn(2), ..., vn(n) do all vanish (including when n = 1)
%     
%     Reference: Golub and Van Loan, "Matrix Computations," 3rd edition,
%     Johns Hopkins University Press, 1996, Chapter 5.
%     

if n==1
    rss = x;
    vn = 1; 
    scal = 0; 
else
    sum = x(2:n)'*x(2:n); 
    
    if sum == 0
       rss = x(1); 
       vn(1) = 1; vn(2:n) = zeros(n-1,1); 
       scal = 0; 
    else
       rss = x(1)^2 + sum;
       rss = sqrt(rss); 
       
       % Determine the first component v1 of the unnormalized Householder vector
       % v = x - phase(x1) * rss * (1 0 0 ... 0 0)^T.
       
       vn(1) = 1; 
       
       if  x(1)<=0
           v1 = x(1) - rss; 
       else
           v1 = -sum / (x(1)+rss); 
       end
       
       %{
       Compute the vector vn and the scalar scal such that vn(1) = 1
       in the Householder transformation identity_matrix - scal * vn * adjoint(vn).
       %}
       
       vn(2:n) = (1/v1)*x(2:n);
       scal = 2*v1^2 / (v1^2+sum);
       
    end
end

