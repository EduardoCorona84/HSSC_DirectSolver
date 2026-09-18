%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%     Test code for HSS1D basic routines: 
%     - compress
%     - invert
%     - transforminv
%     - apply
%     - applyinv
%     
%     We take I_sk (boundary layers) from one box of a uniform binary tree, 
%     and compress A = K[I_sk,I_sk] for some kernel K.  
%     
%     We produce a compressed inverse for A^{-1} in the HSSinv and HSS 1D formats
%     
%     Finally, we test them by using apply and applyinv (for HSS1D and HSSinv
%     structures, respectively) using N random right hand sides and comparing
%     against the corresponding dense applies. 
%     

% 2D Tree indices / points
Nbox = 2;  
lev = TREE.BOX(Nbox).levbox;
I_sk = TREE.BOX(Nbox).I_sk;
X_sk = X(I_sk,:);
k = size(X_sk,1); 
sym = params.sym; 

% Parameters
h = params.h;
lay = params.layers;
n_max = 56; 

% Dense matrix evaluation
A = Kernel_Eval(X_sk,X_sk,params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HSS1D compress (green trick)
tic; 
if sym == 1
    %A_HSS = HSS1D_compress_green_fsym(X_sk,n_max,acc,'Lap',params,1);
    fprintf('\n HSS1D compression fsym, green trick \n')
    A_HSS  = HSS1D_compress_fsym('green',X_sk,n_max,acc,params); 
    fprintf('\n HSS1D compression fsym, brute force \n')
    Ab_HSS = HSS1D_compress_fsym('brute',A,n_max,acc); 
else
    %A_HSS = HSS1D_compress_green_nsym(X_sk,n_max,acc,'Lap',params,1);
    fprintf('\n HSS1D compression nsym, green trick \n')
    A_HSS  = HSS1D_compress_nsym('green',X_sk,n_max,acc,params); 
    fprintf('\n HSS1D compression nsym, brute force \n')
    Ab_HSS = HSS1D_compress_nsym('brute',A,n_max,acc); 
end
T_comp_sk = toc; 

fprintf('\n HSS structure for A (green and brute) \n')
A_HSS(1:9,:)
Ab_HSS(1:9,:)

% Test for N random vector applies
N = 10; 
rand( 'seed',0);

qq_sk = randn(k,N);  

fprintf('\n HSS1D apply \n')
% A*qq
if sym==1
    uHSS_sk  = HSS1D_apply_fsym(A_HSS,qq_sk);
    ubHSS_sk = HSS1D_apply_fsym(Ab_HSS,qq_sk);
else
    uHSS_sk  = HSS1D_apply_nsym(A_HSS,qq_sk);
    ubHSS_sk = HSS1D_apply_nsym(Ab_HSS,qq_sk);
end

% true A*qq
utrue_sk = A*qq_sk; 
    
% Relative Error
Err_Apply_sk  = norm(utrue_sk - uHSS_sk) / norm(utrue_sk); 
Err_Apply_bsk = norm(utrue_sk - ubHSS_sk) / norm(utrue_sk); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HSS1D invert (HSSinv structure) and applyinv
fprintf('\n HSS1D invert and applyinv \n')
if sym==1
    [HSSinv,FTinv] = HSS1D_invert_fsym(A_HSS,acc);
    YHSS_sk = HSS1D_applyinv_fsym(A_HSS,HSSinv,FTinv,qq_sk); 
else
    [HSSinv,FTinv] = HSS1D_invert_nsym(A_HSS,acc);
    YHSS_sk = HSS1D_applyinv_nsym(A_HSS,HSSinv,FTinv,qq_sk); 
end
    
Err_InvApply_sk = norm(YHSS_sk - A\qq_sk) / norm(A\qq_sk); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HSS1D transforminv (from HSSinv to HSS1D structure) and apply
fprintf('\n HSS1D transforminv and apply \n')
if sym==1
    Ainv_HSS = HSS1D_transforminv_fsym(A_HSS,HSSinv,FTinv);
    Y2HSS_sk = HSS1D_apply_fsym(Ainv_HSS,qq_sk); 
else
    Ainv_HSS = HSS1D_transforminv_nsym(A_HSS,HSSinv,FTinv);
    Y2HSS_sk = HSS1D_apply_nsym(Ainv_HSS,qq_sk); 
end

Err_InvApply2_sk = norm(Y2HSS_sk - A\qq_sk) / norm(A\qq_sk); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nMatrix Apply Errors (green and brute):\n')
display(Err_Apply_sk)
display(Err_Apply_bsk)

fprintf('\n Inverse Apply (applyinv HSSinv) Error: \n')
display(Err_InvApply_sk)

fprintf('\n Inverse Apply (apply after transforminv) Error: \n')
display(Err_InvApply2_sk)

