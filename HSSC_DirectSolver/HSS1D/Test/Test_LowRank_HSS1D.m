%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
% 	Test code for HSS1D compression of a low rank matrix. We take T_up from one box
% 	of a HSS2D binary tree (already in low rank form) and compress it into
% 	HSS form. We then test these by applying them to N random
% 	vectors and measuring timings and error. 
%

Test_setup;
rand( 'seed',0);

M = 800; 
k = 50; 

U = randn(M,k); U = U./norm(U); 
V = randn(M,k); V = V./norm(V); 

n_max = 2*k; 

% Conversion to an HSS matrix
tic; 
T_HSS = HSS1D_compressLR_nsym(U,V,n_max,acc);
T_comp_LR = toc; 


% Test for N random vector applies
N = 10; 
 
qq = randn(M,N); 
        
% A*qq
uHSS = HSS1D_apply_nsym(T_HSS,qq,1);
% true A*qq
utrue = U*(V'*qq); 
    
% Relative Error
Err_Apply = norm(utrue - uHSS) / norm(utrue_sk); 
    
fprintf('\n LR to HSS Error \n')
mean(Err_Apply)