function T = LRID(A,acc)
%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
% 	This function writes A symmetric as low rank using the ID. 

[Tid,Jid] = ID(A,acc); k = size(Tid,1);

% U is defined as the subsampled set of columns of A. 
U = A(:,Jid(1:k));

% V' is in turn R, the right interpolation matrix, including the inverse permutation
V(Jid,:) = [eye(k) Tid].'; 

T.U = U; 
T.V = V; 