function [EDvC,ED,EC,Ttr,TDI,TCI,TDA,TCA,MBD,MBC] = Test_HSS2D_CMZ2012_paper(Ng,np)

% Run tests corresponding to the numerical results section of the 
% "An O(N) Direct Solver for Integral Equations on the Plane" 
% by Corona, Martinsson and Zorin.

EDvC = zeros(1,6); ED = EDvC; EC = ED; Ttr = ED; TDI = ED; TCI = ED; 
TDA = ED; TCA = ED; MBD = ED; MBC = ED; 

% Test setup (addpaths, global variables) 
Test_setup; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(1.1) Laplace 2D, acc=1e-10, Translation Invariant kernel

%Parameters
flag_pot = 'SL_L_2D'; kh = 0; sym = 1; acc = 1e-10; lay = 2; n_cut = 400; TI = 1;      
%Generate params struct    
paramsL = HSS_tree_parameters(flag_pot,kh,sym,Ng,np,acc,lay,n_cut,TI); 
paramsL.skel_rule = 'bdry'; paramsL.proxy_rule = 'bdry'; 
paramsL.skel_extra_pts = 0;    
display(paramsL)    

% Test Inverse Compression and Apply
fprintf('\n ------------------------------------------------------------------')
fprintf('\n (5.1) Laplace 2D kernel, TI case, acc=1e-10 \n')
fprintf(' ------------------------------------------------------------------ \n')
[EDvC(1),ED(1),EC(1),Ttr(1),TDI(1),TCI(1),TDA(1),TCA(1),MBD(1),MBC(1)]...
    = Test_HSS2D_InverseCompression(paramsL); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(1.2) Laplace 2D, acc=1e-10, Non Translation Invariant kernel
paramsL.transinv = 0; paramsL.skel_extra_pts = 1; 

% Test Inverse Compression and Apply
fprintf('\n ------------------------------------------------------------------')
fprintf('\n (5.1) Laplace 2D kernel, NTI case, acc=1e-10 \n')
fprintf(' ------------------------------------------------------------------ \n')
[EDvC(2),ED(2),EC(2),Ttr(2),TDI(2),TCI(2),TDA(2),TCA(2),MBD(2),MBC(2)]...
    = Test_HSS2D_InverseCompression(paramsL); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(1.3) Laplace 2D, acc=1e-6, Translation Invariant kernel
paramsL.transinv = 1; paramsL.acc = 1e-6; paramsL.skel_extra_pts = 0; 

% Test Inverse Compression and Apply
fprintf('\n ------------------------------------------------------------------')
fprintf('\n (5.2) Laplace 2D kernel, TI case, acc=1e-6 \n')
fprintf(' ------------------------------------------------------------------ \n')
[EDvC(3),ED(3),EC(3),Ttr(3),TDI(3),TCI(3),TDA(3),TCA(3),MBD(3),MBC(3)]...
    = Test_HSS2D_InverseCompression(paramsL); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(1.4) Laplace 2D, acc=1e-6, Non Translation Invariant kernel
paramsL.transinv = 0; paramsL.skel_extra_pts = 1; 

% Test Inverse Compression and Apply
fprintf('\n ------------------------------------------------------------------')
fprintf('\n (5.2) Laplace 2D kernel, NTI case, acc=1e-6 \n')
fprintf(' ------------------------------------------------------------------ \n')
[EDvC(4),ED(4),EC(4),Ttr(4),TDI(4),TCI(4),TDA(4),TCA(4),MBD(4),MBC(4)]...
    = Test_HSS2D_InverseCompression(paramsL); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(1.5) Helmholtz 2D, acc=1e-10, Translation Invariant kernel
% We generate the parameters string again. This automatically computes the
% Duan-Rokhlin weights for 4th order, for the corresponding k = 8*pi 
% (kappa = 8) 

%Parameters
flag_pot = 'SL_H_2D'; kh = 8*pi; acc = 1e-10; lay = 3; n_cut = 800; TI = 1;        
%Generate params struct    
paramsH = HSS_tree_parameters(flag_pot,kh,sym,Ng,np,acc,lay,n_cut,TI); 
paramsH.skel_rule = 'bdry'; paramsH.proxy_rule = 'rand'; paramsH.lambda = 2; 
paramsH.skel_extra_pts = 1;    
display(paramsH)    

% Test Inverse Compression and Apply
fprintf('\n ------------------------------------------------------------------')
fprintf('\n (5.3) Helmholtz 2D kernel, TI case, acc=1e-10 \n')
fprintf(' ------------------------------------------------------------------ \n')
[EDvC(5),ED(5),EC(5),Ttr(5),TDI(5),TCI(5),TDA(5),TCA(5),MBD(5),MBC(5)]...
    = Test_HSS2D_InverseCompression(paramsH); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(1.6) Lippman-Schwnger 2D, acc=1e-10, Non Translation Invariant kernel
% Again, this is for k=8*pi / kappa = 8. 
paramsH.transinv = 0; 

% Test Inverse Compression and Apply
fprintf('\n ------------------------------------------------------------------')
fprintf('\n (5.4) Lippmann Schwinger (Helmholtz) scattering problem, NTI, acc=1e-10 \n')
fprintf(' ------------------------------------------------------------------ \n')
[EDvC(6),ED(6),EC(6),Ttr(6),TDI(6),TCI(6),TDA(6),TCA(6),MBD(6),MBC(6)]...
    = Test_HSS2D_InverseCompression(paramsH); 

end
