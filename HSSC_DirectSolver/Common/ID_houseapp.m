function [v,scal] = ID_houseapp(n,vn,u,ifrescal,scal)
%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%      Applies the Householder matrix identity_matrix - scal * vn * adjoint(vn)
%      to the vector u, yielding the vector v;
%      
%      scal = 2/(1 + |vn(2)|^2 + ... + |vn(n)|^2) when vn(2), ..., vn(n) don't all vanish;
%      
%      scal = 0 when vn(2), ..., vn(n) do all vanish (including when n = 1).
%      
%      Input:
%      n -- size of vn, u, and v, though the indexing on vn goes from 2 to n
%      vn -- components 2 to n of the Householder vector vn; vn(1) is assumed to be 1 
%      u -- vector to be transformed
%      ifrescal -- set to 1 to recompute scal from vn(2), ..., vn(n);
%                         set to 0 to use scal as input
%      scal -- see the entry for ifrescal in the decription of the input
%      
%      Output:
%      scal -- see the entry for ifrescal in the decription of the input
%      v -- result of applying the Householder matrix to u;
%      
%      Reference: Golub and Van Loan, "Matrix Computations," 3rd edition,
%      Johns Hopkins University Press, 1996, Chapter 5.
%

if n==1
    v = u; 
else
    if ifrescal==1
        %Calculate |vn(2)|^2 + ... + |vn(n)|^2.
        
        sum = vn(2:n)'*vn(2:n); 
        
        if sum == 0
            scal = 0; 
        else
            scal = 2/(1+sum); 
        end
    end
        
    % Calculate fact = scal * vn' * u.
    fact = scal*(vn'*u); 
    % Subtract fact*vn from u, yielding v
    v = u - fact*vn; 

end

