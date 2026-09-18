function Y = BHSS1D_apply_fsym(K11,K12,K21,K22,B_dn,B_up,I_dn,I_up,X) 
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        Y = BHSS1D_apply_fsym(K11,K12,K21,K22,B_dn,B_up,I_dn,I_up,X)
%
%    DESCRIPTION:
%        This function performs a fast matrix-vector multiply Y=AX, where A is a submatrix of a 2x2
%        block matrix with full-symmetric HSS1D elements; it is specified by parameters
%        {K11,K12,K21,K22,B_dn,B_up,I_dn,I_up}. X is a given input vector and Y is the solution
%        vector.
%
%    INPUT:
%        A denotes the submatrix of the 2x2 HSS1D block matrix
%            K = [K11 K12 ; K21 K22]
%        I_dn / I_up {i,j} are boolean index vectors that specify the submatrix A as follows:
%            A = [K11(I_dn{1,1},I_up{1,1}) K12(I_dn{1,2},I_up{1,2}) ;
%                 K21(I_dn{2,1},I_up{2,1}) K22(I_dn{2,2},I_up{2,2})]
%        B_dn / B_up {i,j} are boolean index vectors that specify the boxes in the subtree taken
%            from Kij. They can to be precomputed as follows:
%            [B_dn{i,j},B_up{i,j}] = HSS1D_apply_submatrix_boxes(Kij,I_dn{i,j},I_up{i,j});
%        X = [ X1 ; X2 ] is the vector to be left multiplied by A.
%
%    OUTPUT:
%        Y = AX is the vector with
%            [Y1 ; Y2] = [ A11*X1 + A12*X2 ; A21*X1 + A22*X2 ];
%


q = size(X,2); 

% submatrix block sizes
k1=sum(I_dn{1,1}); r1=sum(I_up{1,1}); n1 = k1 + r1;
k2=sum(I_dn{2,2}); r2=sum(I_up{2,2}); n2 = k2 + r2;

% Separate X1 and X2
X1 = X(1:r1,:); 
X2 = X(r1+1:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vec(I_up,:) = X1, zero elsewhere
Vec1 = zeros(n1,q); 
Vec1(I_up{1,1},:) = X1;
% Y11 = K11*Vec = A11*X1
Y11 = HSS1D_apply_submatrix_fsym(K11,B_dn{1,1},B_up{1,1},I_dn{1,1},I_up{1,1},Vec1,1);

% Vec(I_up,:) = X2, zero elsewhere
Vec1 = zeros(k1+r2,q); 
Vec1(I_up{1,2},:) = X2;
% Y12 = K12*Vec = A12*X2
Y12 = HSS1D_apply_submatrix_fsym(K12,B_dn{1,2},B_up{1,2},I_dn{1,2},I_up{1,2},Vec1,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vec(I_up,:) = X1, zero elsewhere
Vec2 = zeros(k2+r1,q);
Vec2(I_up{2,1},:) = X1;
% Y21 = K21*Vec = A21*X1
Y21 = HSS1D_apply_submatrix_fsym(K21,B_dn{2,1},B_up{2,1},I_dn{2,1},I_up{2,1},Vec2,1);

% Vec(I_up,:) = X2, zero elsewhere
Vec2 = zeros(n2,q); 
Vec2(I_up{2,2},:) = X2;
% Y22 = K22*Vec = A22*X2
Y22 = HSS1D_apply_submatrix_fsym(K22,B_dn{2,2},B_up{2,2},I_dn{2,2},I_up{2,2},Vec2,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Y

Y = zeros(k1+k2,q);

Y(1:k1,:)     = Y11 + Y12;
Y(k1+1:end,:) = Y21 + Y22;