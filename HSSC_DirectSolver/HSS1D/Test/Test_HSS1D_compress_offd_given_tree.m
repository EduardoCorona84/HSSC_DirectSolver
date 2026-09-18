%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%   Test code for HSS1D compression given a tree. We take I_sk from one box
%   of a uniform binary tree, sort it and compress A_sk = K[I_sk,I_sk]. 
%   We use the resulting tree structure to test compression given a tree,
%   and then test these by applying them to N random
%   vectors and measuring timings and error. 
%

% 2D Tree indices / points
Nbox = 3; tic; toc; 
X = params.X_source;
lev = TREE.BOX(Nbox).levbox;
I_sk = TREE.BOX(Nbox).I_sk;
X_sk = X(I_sk,:);
k = size(X_sk,1);  

% parameters
h = params.h;
lay = params.layers;
cB = TREE.BOX(Nbox).cent;
n_max = 50; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HSS1D matrix compression (green offd vs offd given tree)

if sym == 1
    tic; 
    %A_off = HSS1D_compress_green_fsym_offd(X_sk,n_max,acc/100,'Lap',params,1);
    A_off  = HSS1D_compress_fsym('green',X_sk,n_max,acc/100,params,'offd'); 
    T_comp_sk = toc; 
    
    tic; 
    %A_off_gt = HSS1D_compress_fsym_offd_giventree(X_sk,A_off,acc/100,'Lap',params,1);
    A_off_gt  = HSS1D_compress_fsym('green',X_sk,A_off,acc/100,params,'offd'); 
    T_comp_gt = toc; 
else
    tic; 
    %A_off = HSS1D_compress_green_nsym_offd(X_sk,n_max,acc/100,'Lap',params,1);
    A_off  = HSS1D_compress_nsym('green',X_sk,n_max,acc/100,params,'offd');
    T_comp_sk = toc; 
    
    tic; 
    %A_off_gt = HSS1D_compress_nsym_offd_giventree(X_sk,A_off,acc/100,'Lap',params,1); 
    A_off_gt  = HSS1D_compress_nsym('green',X_sk,A_off,acc/100,params,'offd'); 
    T_comp_gt = toc;
end

% Dense evaluation of K(X,X)
K = Kernel_Eval(X_sk,X_sk,params); 
fprintf('\n Is K symmetric? \n')
norm(K - K.')

% Dense evaluation of off-diagonal 2 x 2 blocks
M = A_off{7,2}; 
M2 = length(K) - M; 
K_off = K; 
K_off(1:M,1:M) = zeros(M,M); 
K_off(M+1:end,M+1:end) = zeros(M2,M2); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compress A and delete necessary blocks
if sym == 1
    %A = HSS1D_compress_green_fsym(X_sk,n_max,acc/100,'Lap',params,1);
    A  = HSS1D_compress_fsym('green',X_sk,n_max,acc/100,params); 
else
    %A = HSS1D_compress_green_nsym(X_sk,n_max,acc/100,'Lap',params,1);
    A  = HSS1D_compress_nsym('green',X_sk,n_max,acc/100,params); 
end

nbox = size(A,2); 
for ibox = 2:nbox
    if ( (A{BOX.C1,ibox}>0) && (A{BOX.C2,ibox}>0) )
        ison1 = A{BOX.C1,ibox}; ison2 = A{BOX.C2,ibox}; 
    A{BOX.M_SIB,ison1} = zeros(size(A{BOX.M_SIB,ison1})); A{BOX.M_SIB,ison2} = A{BOX.M_SIB,ison1}.'; 
    else
       A{BOX.M_SELF,ibox} = zeros(size(A{BOX.M_SELF,ibox}));  
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Test for N random vector applies
N = 10; 
rand( 'seed',0);
qq_sk = randn(k,N); qq_sk = qq_sk./norm(qq_sk); 
 
% A*qq
if sym == 1
    uHSS_sk = HSS1D_apply_fsym(A,qq_sk);
    uHSS_gt = HSS1D_apply_fsym(A_off_gt,qq_sk);
else
    uHSS_sk = HSS1D_apply_nsym(A,qq_sk);
    uHSS_gt = HSS1D_apply_nsym(A_off_gt,qq_sk);
end

u = K_off*qq_sk; 
    
% Relative Error
Err_Apply = norm(uHSS_gt - uHSS_sk) / norm(uHSS_sk); 
Err_gt    = norm(uHSS_gt - u) / norm(u); 

fprintf('\n Off-d Error (HSS vs Dense) \n')
display(Err_gt)

fprintf('\n Off-d Error (off-d vs compress + 0 off-d blocks) \n')
display(Err_Apply)
