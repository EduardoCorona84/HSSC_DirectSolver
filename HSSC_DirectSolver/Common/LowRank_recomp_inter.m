function [U,V] = LR_recompInter(U_cell,V_cell,acc)
%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%	  Function that, given cells: 
%	  U_cell = {U1,U2,...,Up} 
%	  V_cell = {V1,V2,...,Vp}
%	  Computes a Low Rank decomposition of the matrix: 
%	  
%	  A = [U1*V1' U2*V2' ... Up*Vp']
%	  
%	  It outputs U and V such that A = U*V'. 
%	  
%	  Inputs: 
%	  U_cell (cell array L x 1) U matrices
%	  V_cell (cell array L x 1) V matrices
%	  acc    (int) accuracy for IDs
%	  
%	  Outputs: 
%	  U      (array m x k) 
%	  V      (array n x k) 
%	  
%	  This is useful in compressing the interpolation operators, which are
%	  computed through low rank translation operators from well-separated pieces 
%	  of the interface to the boundary of the box via proxy rectangles. 
%	  
%	  We note that here, the U matrices can be compressed together, but the V
%	  matrices cannot (since the pieces need not be the same size). 
%

L = length(U_cell);
 
U_mat = U_cell{1}'; 

for i=2:L
   U_mat = [U_mat ; U_cell{i}']; 
end

[T_dn J_dn] = ID(U_mat,acc); k_dn = size(T_dn,1); 

% U is defined as the Left interpolation matrix (including the inverse
% permutation)
U(J_dn,:) = [eye(k_dn) ; T_dn']; 

% V' is in turn S*R, where S is A(J_dn(1:k_dn),J_up(1:k_up)) and R is the right 
% interpolation matrix, again including the inverse permutation

S = U_cell{1}(J_dn(1:k_dn),:)*V_cell{1}'; 

for i = 2:L
    S = [S U_cell{i}(J_dn(1:k_dn),:)*V_cell{i}']; 
end

V = S'; 
    