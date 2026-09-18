function [Z,k] = Rand_AT_sample(m,n,matvecA,matvecAt,dist,params,k0,C)
%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%   This function implements the randomized version of pivoted QR. 
%   A partial double Gram-Schmidt is used to find the numerical rank k
%   and the orthogonal matrix Q. The goal is to find an orthonormal basis for
%   the range of A. 
%   
%   
%   suggested value for tol
%   tol = eps*norm(A,'fro');
%   If a value of k is not provided a priori, we compute the numerical rank 
%   of A using the double GS
%   
%   Programmer: Eduardo Corona
%   Numerical Linear Algebra with Probability, Mark Tygert 2010
%


if nargin == 8
    tol = params; 
    [k,Y] = QR_findrank(m,n,matvecAt,dist,k0,C,tol);
    %[k,Y] = ID_findrank(m,n,matvecAt,dist,k0,tol);
    Z = Y.';
else
    l=params; k=l;
    
    % Random Matrix(Gaussian or Uniform depending on dist)
    if strcmp(dist,'unif')
        G = rand(m,l);
    else
        G = randn(m,l);
    end

    Y = matvecAt(G); 
    Z = Y.';
end


