function [U,V] = LR_recomp(U_cell,V_cell,acc)
%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%		Function that, given cells: 
%		U_cell = {U1,U2,...,Up} 
%		V_cell = {V1,V2,...,Vp}
%		Computes a Low Rank decomposition of the matrix: 
%		
%		A = U1*V1' + U2*V2' + ... + Up*Vp'
%		
%		It outputs U and V such that A = U*V'. 
%		
%		Inputs: 
%		U_cell (cell array L x 1) U matrices
%		V_cell (cell array L x 1) V matrices
%		acc    (int) accuracy for IDs
%		
%		Outputs: 
%		U      (array m x k) 
%		V      (array n x k) 
%		


L = length(V_cell);

V_mat = V_cell{1}.'; 

for i=2:L
   V_mat = [V_mat ; V_cell{i}.' ]; 
end

[T_up, J_up] = ID(V_mat,acc); k_up = size(T_up,1); 

U = U_cell{1}*(V_cell{1}(J_up(1:k_up),:)).';
% U is defined as the subsampled set of columns of A. 
for i = 2:L
    Utmp = U_cell{i}; Vtmp = V_cell{i}.'; 
    U = U + Utmp*Vtmp(:,J_up(1:k_up)); 
end

% V' is in turn R, the right interpolation matrix, including the inverse permutation
V(J_up,:) = [ eye(k_up) T_up].'; 

clear Utmp Vtmp V_mat; 
    