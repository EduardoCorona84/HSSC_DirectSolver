%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%   Test code for HSS1D sum and recompression. We take I_sk from one box
%   of a HSS2D binary tree, sort it and compress 
%   
%   A = K1[I_sk,I_sk] and B = K2[I_sk,I_sk] 
%   
%   for two different kernel definitions. 
%   
%   We then compare the applies of the sum (obtained with HSS1D_sum) and the
%   sum of the two HSS applies for A and B for N random vectors, measuring 
%   timings and error. 
%   
%   Finally, we look at the HSS1D src and skeleton set sizes before and after
%   recompression. 
%


% 2D Tree indices / points
Nbox = 7; lev = 2; 
X = params.X_source;
I_sk = TREE.BOX(Nbox).I_sk;
X_sk = X(I_sk,:);
k = size(X_sk,1);
cB = TREE.BOX(Nbox).cent;

% parameters
h = params.h;
lay = params.layers;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A matrix compress
if params.sym == 0
    tic; 
    %A_HSS = HSS1D_compress_green_nsym(X_sk,56,acc,'Lap',params,1);
    A_HSS = HSS1D_compress_nsym('green',X_sk,56,acc,params); 
    T_compA = toc; 
else
    tic; 
    %A_HSS = HSS1D_compress_green_fsym(X_sk,56,acc,'Lap',params,1);
    A_HSS = HSS1D_compress_fsym('green',X_sk,56,acc,params); 
    T_compA = toc; 
end

paramsB = params; paramsB.h = 100*h; 

% B matrix compress
if params.sym == 0
    tic; 
    %B_HSS = HSS1D_compress_green_nsym(X_sk,56,acc,'Lap',paramsB,1);
    B_HSS = HSS1D_compress_nsym('green',X_sk,56,acc,paramsB); 
    T_compB = toc; 
else
    tic; 
    %B_HSS = HSS1D_compress_green_fsym(X_sk,56,acc,'Lap',paramsB,1);
    B_HSS = HSS1D_compress_fsym('green',X_sk,56,acc,paramsB); 
    T_compB = toc; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C = A + B (sum and recompress)
if params.sym == 0
    tic; 
    C_HSS = HSS1D_sum_nsym(A_HSS,B_HSS,acc);
    T_sum = toc; 

    CR_HSS = HSS1D_recompress_nsym(C_HSS,acc);
else
    tic; 
    C_HSS = HSS1D_sum_fsym(A_HSS,B_HSS,acc);
    T_sum = toc; 
    
    CR_HSS = HSS1D_recompress_fsym(C_HSS,acc);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test for N random vector applies
N = 10; 
rand( 'seed',0);

qq = randn(k,N); 
   
if params.sym == 0
    % A*qq
    uA = HSS1D_apply_nsym(A_HSS,qq);
    %B*qq
    uB = HSS1D_apply_nsym(B_HSS,qq);
else
    % A*qq
    uA = HSS1D_apply_fsym(A_HSS,qq);
    %B*qq
    uB = HSS1D_apply_fsym(B_HSS,qq);
end
    
% A*qq + B*qq
u_sum = uA+uB; 
    
if params.sym == 0
    % C*qq
    uC = HSS1D_apply_nsym(C_HSS,qq);
    uCR = HSS1D_apply_nsym(CR_HSS,qq);
else
    % C*qq
    uC = HSS1D_apply_fsym(C_HSS,qq);
    uCR = HSS1D_apply_fsym(CR_HSS,qq);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Measure Relative Error
Err_sum = norm(uC - u_sum) / norm(u_sum);
Err_rec = norm(uC - uCR) / norm(u_sum);

fprintf('\n\n Sum Application Error \n')
display(Err_sum)

fprintf('\n Recompress Error \n')
display(Err_rec)
    
fprintf('\n Recompress Effect on Skeleton Sizes \n')
fprintf('\n Before: \n')
C_HSS(7:9,:)
fprintf('\n After: \n')
CR_HSS(7:9,:)