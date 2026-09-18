%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%   Test code for HSS1D split into two diagonal sub-blocks K11 and K22. . 
%   We take I_sk from one box of a HSS2D binary tree, sort it and compress 
%   A_HSS = K[I_sk,I_sk]. 
%   
%   Given a boolean index B, HSS1D split produces HSS forms for K[B,B] and 
%   K[~B,~B]. 
%   
%   We then test the split routine for 
%   (1) half and half
%   (2) 1/4 vs 3/4 
%   (3) even vs odd
%   (4) randomly picked with zeroed out beginning and ending
%   
%   Finally, we test these by applying them to N random
%   vectors and measuring timings and error. 
%

% 2D Tree indices / points
Nbox = 7; 
X = params.X_source;
lev = TREE.BOX(Nbox).levbox;
I_sk = TREE.BOX(Nbox).I_sk;
X_sk = X(I_sk,:);
k = size(X_sk,1); 
cB = TREE.BOX(Nbox).cent;

% Parameters
h = params.h;
lay = params.layers;
n_max = 56; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HSS matrices compress
if sym == 1
tic; 
%A_HSS = HSS1D_compress_green_fsym(X_sk,n_max,acc,'Lap',params,1);
A_HSS  = HSS1D_compress_fsym('green',X_sk,n_max,acc,params);
T_comp_sk = toc;
else
    tic; 
%A_HSS = HSS1D_compress_green_nsym(X_sk,n_max,acc,'Lap',params,1);
A_HSS  = HSS1D_compress_nsym('green',X_sk,n_max,acc,params);
T_comp_sk = toc;
end

% Different boolean vectors to split HSS matrix

I1 = (1:k)<(k/2+1); N1 = sum(I1); 
I2 = (1:k)<(k/4+1); N2 = sum(I2); 
I3 = mod((1:k),2)==0 ; N3 = sum(I3); 

% Randomly picked points
Irand = rand(1,k); 
pr = 0.8; 
for i=1:k
   if Irand(i)>(1-pr)
       Irand(i) = 1; 
   else
       Irand(i) = 0; 
   end
end
% Blank sections
Irand(1:39) = zeros(1,39); 
Irand(k-37:end) = zeros(1,38); 
Irand = Irand==1 ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HSS1D split for 4 cases

if sym == 1
% 1/2 - 1/2
tic;
[KL1,KR1] = HSS1D_split_fsym(A_HSS,I1,acc); 
T_split1 = toc; 
% 1/4 - 3/4

tic;
[KL2,KR2] = HSS1D_split_fsym(A_HSS,I2,acc); 
T_split2 = toc; 

% even points
tic;
[KL3,KR3] = HSS1D_split_fsym(A_HSS,I3,acc); 
T_split3 = toc; 

% at random plus zeros
tic;
[KL4,KR4] = HSS1D_split_fsym(A_HSS,Irand,acc); 
T_split4 = toc;
else
    % 1/2 - 1/2
tic;
[KL1,KR1] = HSS1D_split_nsym(A_HSS,I1,I1,acc); 
T_split1 = toc; 
% 1/4 - 3/4

tic;
[KL2,KR2] = HSS1D_split_nsym(A_HSS,I2,I2,acc); 
T_split2 = toc; 

% even points
tic;
[KL3,KR3] = HSS1D_split_nsym(A_HSS,I3,I3,acc); 
T_split3 = toc; 

% at random plus zeros
tic;
[KL4,KR4] = HSS1D_split_nsym(A_HSS,Irand,Irand,acc); 
T_split4 = toc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test for N random vector applies
N = 10; 
rand( 'seed',0);

qq_1 = randn(k,N); 
qq_2 = randn(k,N);
qq_3 = randn(k,N); 
qq_4 = randn(k,N); 

qq_even = zeros(k,N); 
qq_even(I3,:) = qq_3(I3,:); 

qq_r = zeros(k,N); 
qq_r(Irand,:) = qq_4(Irand,:); 

if sym == 1    
    % A*qq
    U1 = HSS1D_apply_fsym(A_HSS,[qq_1(1:N1,:) ; zeros(k-N1,N)]);
    U2 = HSS1D_apply_fsym(A_HSS,[qq_2(1:N2,:) ; zeros(k-N2,N)]);
    U3 = HSS1D_apply_fsym(A_HSS,qq_even);
    U4 = HSS1D_apply_fsym(A_HSS,qq_r);
    
    UL1 = HSS1D_apply_fsym(KL1,qq_1(1:N1,:));
    UL2 = HSS1D_apply_fsym(KL2,qq_2(1:N2,:));
    UL3 = HSS1D_apply_fsym(KL3,qq_3(I3,:));
    UL4 = HSS1D_apply_fsym(KL4,qq_4(Irand,:)); 
else
    % A*qq
    U1 = HSS1D_apply_nsym(A_HSS,[qq_1(1:N1,:) ; zeros(k-N1,N)]);
    U2 = HSS1D_apply_nsym(A_HSS,[qq_2(1:N2,:) ; zeros(k-N2,N)]);
    U3 = HSS1D_apply_nsym(A_HSS,qq_even);
    U4 = HSS1D_apply_nsym(A_HSS,qq_r);
    
    UL1 = HSS1D_apply_nsym(KL1,qq_1(1:N1,:));
    UL2 = HSS1D_apply_nsym(KL2,qq_2(1:N2,:));
    UL3 = HSS1D_apply_nsym(KL3,qq_3(I3,:));
    UL4 = HSS1D_apply_nsym(KL4,qq_4(Irand,:)); 
end
    
    % Relative Error
    Err_Apply_L1 = norm(U1(I1,:) - UL1) / norm(U1(I1,:)); 
    Err_Apply_L2 = norm(U2(I2,:) - UL2) / norm(U2(I2,:)); 
    Err_Apply_L3 = norm(U3(I3,:) - UL3) / norm(U3(I3,:)); 
    Err_Apply_L4 = norm(U4(Irand,:) - UL4) / norm(U4(Irand,:)); 


fprintf('Error Apply 1 (half and half)'); 
mean(Err_Apply_L1)

fprintf('Error Apply 2 (1/4 and 3/4)'); 
mean(Err_Apply_L2)

fprintf('Error Apply 3 (even and odd)'); 
mean(Err_Apply_L3)

fprintf('Error Apply 4 (random and ends off)'); 
mean(Err_Apply_L4)