% Run All HSS2D tests

Test_setup;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters
flag_pot = 'SL_L_2D'; kh = 0; sym = 1; Ng = 16; np = 8; acc = 1e-10; lay = 2; n_cut = 400; TI = 1;      
%Generate params struct    
params = HSS_tree_parameters(flag_pot,kh,sym,Ng,np,acc,lay,n_cut,TI); 
params.skel_rule = 'bdry'; params.dim = 2; params.skel_extra_pts = 0;    
display(params)    

% Test binary tree
fprintf('\n ------------------------------------------------------------------')
fprintf('\n Test_HSS2D_bintree\n')
fprintf(' ------------------------------------------------------------------ \n')
Test_HSS2D_bintree;

% Test Interpolation Operator routines
%fprintf('\n\n Test_interpolation_operator\n\n')
%Test_interpolation_operator;

% Test Matrix Compression and Apply
fprintf('\n ------------------------------------------------------------------')
fprintf('\n Test Matrix Compress and Apply \n')
fprintf(' ------------------------------------------------------------------ \n')
Test_HSS2D_Matrix_Apply(params,INFO_TREE);

% Test Inverse Compression and Apply
fprintf('\n ------------------------------------------------------------------')
fprintf('\n Test Inverse Compression \n')
fprintf(' ------------------------------------------------------------------ \n')
Test_HSS2D_InverseCompression(params,INFO_TREE); 