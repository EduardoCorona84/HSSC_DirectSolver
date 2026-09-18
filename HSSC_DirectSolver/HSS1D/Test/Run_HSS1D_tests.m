%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%	 This Test routine runs all relevant HSS1D arithmetic tests
%

st = clock;
Test_setup;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters
flag_pot = 'SL_L_2D'; kh = 0; sym = 0; Ng = 7; np = 16; acc = 1e-10; n_cut = 432; lay = 2; TI = 0; 
%Generate params struct 
params = HSS_tree_parameters(flag_pot,kh,sym,Ng,np,acc,lay,n_cut,TI);
params.dim = 2; 
params.proxy = 1; 

display(params)

% Build 2D TREE structure
fprintf('\n ------------------------------------------------------------------')
fprintf('\n Build 2D binary tree on uniform grid of [-1,1]^2 \n')
fprintf(' ------------------------------------------------------------------ \n')
TREE = Unifbintree(params);
X = params.X_source;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OMNI routines: compress, invert and applies
fprintf('\n ------------------------------------------------------------------')
fprintf('\n Test HSS1D compression, invert, transforminv, apply and applyinv \n')
fprintf(' ------------------------------------------------------------------ \n')
Test_HSS1D_compress; 

% Off-diagonal compression (green_offd and offd_given_tree)
fprintf('\n ------------------------------------------------------------------')
fprintf('\n Test HSS1D compress offd given tree \n')
fprintf(' ------------------------------------------------------------------ \n')
Test_HSS1D_compress_offd_given_tree; 

% Sum of HSS matrices (add test for Recompress)
fprintf('\n ------------------------------------------------------------------')
fprintf('\n Test HSS1D sum and recompress \n')
fprintf(' ------------------------------------------------------------------ \n')
Test_HSS1D_SUM;

% Split test for multiple cases
fprintf('\n ------------------------------------------------------------------')
fprintf('\n Test HSS1D split \n')
fprintf(' ------------------------------------------------------------------ \n')
Test_HSS1D_split;

% Series of splits, merges, offd, sums imitating Build_FInv
fprintf('\n ------------------------------------------------------------------')
fprintf('\n Test HSS1D compress offd, split, merge and sum \n')
fprintf(' ------------------------------------------------------------------ \n')
Test_HSS1D_merge_and_sum(TREE,params); 

% Low Rank to HSS1D
fprintf('\n ------------------------------------------------------------------')
fprintf('\n Test HSS1D compress of Low Rank operator \n')
fprintf(' ------------------------------------------------------------------ \n')
Test_LowRank_HSS1D

fprintf(' ------------------------------------------------------------------ \n')
fprintf(' Total Time: ')
disp(etime(clock, st))
