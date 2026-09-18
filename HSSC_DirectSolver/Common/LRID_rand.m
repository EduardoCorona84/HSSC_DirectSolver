function T = LRID_rand(m,n,A,At,dist,params)
%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%	 This function writes A as low rank using the ID. 

if strcmp(params.form,'HSS')
    matvecA = @(x) OMNI_Apply_HSS_fsym(A,x); 
    matvecAt = @(x) OMNI_Apply_HSS_fsym(At,x);

else
    %Else, use a matvec for block HSS. 
    matvecA = A; 
    matvecAt = At;
end

par = params.par; 
if par<1
    k0 = params.q;
    C  = params.C; 
    [Ti J] = ID_rand(m,n,matvecA,matvecAt,dist,par,k0,C);
else
    [Ti J] = ID_rand(m,n,matvecA,matvecAt,dist,par);
end
    
k = size(Ti,1);

I = eye(n); 
AJ = matvecA(I(:,J(1:k)));

% U is defined as the subsampled set of columns of A. 
U = AJ; 
% V' is in turn R, the right interpolation matrix, including the inverse permutation
V(J,:) = [eye(k) Ti].'; 

T.U = U; 
T.V = V; 