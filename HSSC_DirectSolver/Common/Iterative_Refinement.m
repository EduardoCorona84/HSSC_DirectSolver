function X = Iterative_Refinement(Ainv,A,B,acc,maxit)
%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
% 	Routine to perform adaptive refinement for the system AX=B and tolerance
% 	acc. A is an N x N invertible matrix and B is N x q. 
%	
%	Inputs: 
%	
%	Ainv (function handle) - Approximate Inverse
%	A    (function handle) - Exact or Approximate A matvec
%	B    (N x q array)     - right-hand-side 
%	acc  (double)          - accuracy
%	maxit (int)            - maximum iterations
%

X = Ainv(B); 
Res = B - A(X); 
Err = norm(Res);
iter = 1; 

display(iter)
display(Err)

while Err > acc & iter < maxit
    X = X + Ainv(Res); 
    Res = B - A(X); 
    Err = norm(Res); 
    iter = iter + 1; 
    display(iter)
    display(Err)
end

