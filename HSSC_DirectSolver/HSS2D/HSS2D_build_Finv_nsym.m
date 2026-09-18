function [FI] = HSS2D_build_Finv_nsym(E1, E2, lev, X, J, k, I1_sk, I2_sk, params, opts)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%     FUNCTION CALL:
%         HSS2D_build_Finv_nsym(E1,E2,lev,X,J,k,I1_sk,I2_sk,params,opts)
%
%     DESCRIPTION:
%         This routine returns the four matrices needed to apply F^-1 by the Schur complement
%         formulae. For a non-leaf box, F is a linear operator defined on the source set.
%         Physically, it maps charge distributions to fields on this set, adding contributions
%         from local operators E1,E2 (of the children) and sibling interactions:
%
%                         F = [E1 K12 ; K21 E2]
%
%         where K12,K21 are kernel evaluations. No assumptions are made on the symmetry of the
%         associated kernel (nsym). This function is an implementation of Algorithm 4 in the
%         paper.
%
%         The speed of the algorithm depends on how E and F^-1 are computed/stored, which is
%         specified by the argument opts, or rather user input params.Fopts of parent function
%         HSS2D_Inverse-Compression_nsym (where HSS2D_build_Finv_nsym is called).
%
%         For the O(N) algorithm, we reorder the rows and columns by skeleton and residual sets,
%         partitioning F into the four blocks F_sk, F_r2s, F_s2r, F_rs
%         (F^{sk}, F{s<-r}, F{r<-s}, F{rs} in the paper).
%
%                         F = [F_sk F_r2s ; F_s2r F_rs]
%                         F_* = [E1_* K12_* ; K21_* E2_*]
%
%         Performing block Gaussian elimination we find F^-1 can be represented by the four
%         matrices F_rs^-1, F_r2s, F_s2r, and S^-1 where S is the Schur complement matrix
%
%                         S = F_sk - F_r2s F_rs^-1 F_s2r
%
%         The computation of these matrices is accelerated by HSS1D and low-rank operations.
%         Because the skeleton sets are fixed-width around the boundary, F_sk and F_rs can be
%         interpreted as HSS1D. Because the skeleton and residual sets are close only at a few
%         points, F_s2r and F_r2s are lowrank.
%
%     INPUT:
%         E1,E2           <kxk>           E matrices of the children, stored densely or as HSS1D.
%         lev             <int>           Box level in the tree.
%         X               <2xn float>     X(I_src,:) = [X(I1_sk); X(I2_sk)] Merge of
%                                                 skeleton points of children (source set).
%         J               <nx1 int>       = [I_sk I_rs] Indices of skeleton and residual
%                                                 points.
%         k               <1x2 int>       = [k1 k2] Size of children skeleton sets
%         I1_sk,I2_sk     <kx1 int>       Indices of children skeleton sets.
%         params          <struct>        (see HSS_tree_parameters) Solve parameters.
%
%         opts    <string>    {'dense', 'block', 'd2HSS', 'HSS'} Specifies how to compute the
%                                 blocks of F^-1. It is determined by the lev, params.lev_HSS,
%                                 and params.Fopts in the parent function
%                                 HSS2D_invert_nsym.
%
%                             opts = Fopts
%                             if (lev > lev_HSS) : opts = 'dense'
%                             if (lev = lev_HSS and Fopts = 'HSS') : opts = 'd2HSS'
%
%                             - 'dense'   F^-1 is computed as a single dense matrix. E matrices
%                                             are assumed dense.
%                             - 'block'   The blocks of FI.* are computed densely. This
%                                             corresponds to the left side of Algorithm 4.
%                             - 'd2HSS'   The blocks of F^-1 are HSS1D or low-rank. E matrices
%                                             are assumed dense.
%                             - 'HSS'     The blocks of F^-1 are HSS1D or low-rank. E matrices
%                                             are HSS1D. This corresponds to the right side of
%                                             Algorithm 4.
%                             - 'HSStest' Same as 'HSS' but it computes errors of each operation
%                                             for debugging purposes.
%
%     OUTPUT:
%         FI either stored densely or as a struct of blocks (either dense or compressed).
%             <nxn float>     Dense matrix of F^-1
%             <struct>        .rs_Inv     HSS1D of residual matrix F_rs^-1
%                             .rs2sk      Low-rank decomposition of F_r2s by randomized ID
%                             .sk2rs      Low-rank decomposition of F_s2r by randomized ID
%                             .SrsInv     Inverse of the Schur complement matrix S.
%                             .J          <mx1 int> = J
%




% If opts == 'dense', then FI is built as a single, dense matrix as in the
% 1D code

if strcmp(opts,'dense')
    
    X1 = X(I1_sk,:); X2 = X(I2_sk,:);
        
    K12 = Kernel_Eval(X1,X2,params); 
    K21 = Kernel_Eval(X2,X1,params); 
    
    F = [E1 K12 ; K21 E2]; 
    FI = inv(F); 
    
elseif strcmp(opts,'HSS')
    
    X1 = X(I1_sk,:); X2 = X(I2_sk,:);
    
    % If opts == 'HSS', then FI is built using the Schur complement
    % formulae. It is assumed that E1 and E2 are HSS matrices. 
    % The diagonal blocks of F are computed as E_dg + K_offd using HSS1DSUM. 
    
    acc = params.acc;
    JE1 = params.JE1; JE2 = params.JE2; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % (I) Computing F_sk and F_rs in HSS form: 
    X1_s = X1(JE1,:); X2_s = X2(JE2,:); 
    
    % (0) SORT: we sort points and separate skeleton and residual sets.
    Bl = false(size(J')); Bl(J(1:k)) = true; 
    
    m1 = size(X1,1); 
    B1 = Bl(1:m1); B2 = Bl(m1+1:end); 
    B1 = B1(JE1); B2 = B2(JE2); 
        
    % (1) SPLIT: we split the E HSS matrices into skeleton and residual sets
    [E1_sk,E1_rs] = HSS1D_split_nsym(E1,B1,B1,acc); 
    [E2_sk,E2_rs] = HSS1D_split_nsym(E2,B2,B2,acc); 
    
    % Consolidate points
    Xsk  = [X1_s(B1,:) ; X2_s(B2,:)]; 
    Xrs  = [X1_s(~B1,:) ; X2_s(~B2,:)];  
    
    % (2) MERGE
    E_sk              = HSS1D_merge_nsym(E1_sk,E2_sk,Xsk,lev,acc,'concatenate');
    [E_rs,Xrs_s,Jrs]  = HSS1D_merge_nsym(E1_rs,E2_rs,Xrs,lev,acc,'align');
    
    Brs = [true(1,sum(~B1)) false(1,sum(~B2))]; 
    Brs = Brs(Jrs); 
    
    % (3) OFFD and SUM
    %Ksk_off  = HSS1D_compress_nsym_offd_giventree(Xsk,E_sk,acc,'lap',params,1);
    Ksk_off  = HSS1D_compress_nsym('green',Xsk,E_sk,acc,params,'offd'); 
    Krs_off  = HSS1D_compress_nsym_offd_leaf(Xrs_s,E_rs,Brs,acc,params); 
    
    F_sk  = HSS1D_sum_nsym(E_sk,Ksk_off,acc); 
    F_rs  = HSS1D_sum_nsym(E_rs,Krs_off,acc);
    F_rs = HSS1D_recompress_nsym(F_rs,acc);  
    
    % Invert F_rs using HSS1D inversion
    [Frs_HSSinv,Frs_FTinv] = HSS1D_invert_nsym(F_rs,acc);
    FI.rs_Inv = HSS1D_transforminv_nsym(F_rs,Frs_HSSinv,Frs_FTinv);    
    
    J1 = JE1; J2 = JE2; 
    JFr = [J1(~B1) (J2(~B2) + length(B1))];
    FI.J = [J1(B1) (J2(B2) + length(B1)) JFr(Jrs)];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (II) Off-diagonal blocks of F as low rank using ID_rand
    
    [T_r2s,T_s2r] = LOCAL_LR_Tr2s(E1,E2,B1,B2,X1_s,X2_s,lev,acc,params);
    
    T_r2s.V = T_r2s.V(Jrs,:); T_s2r.U = T_s2r.U(Jrs,:); 
    
    FI.rs2sk = T_r2s; 
    FI.sk2rs = T_s2r;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (III) Schur Complement Matrix S_rs
    
    U2 = HSS1D_apply_nsym(FI.rs_Inv,T_s2r.U); 
    U = -LR_apply(T_r2s,U2); 
    V = T_s2r.V; 
    
    M = HSS1D_compressLR_nsym(U,V,F_sk,acc);
    S_rs = HSS1D_sum_nsym(F_sk,M,acc); 
    S_rs = HSS1D_recompress_nsym(S_rs,acc);
 
    [S_HSSinv,S_FTinv] = HSS1D_invert_nsym(S_rs,acc);
    FI.SrsInv = HSS1D_transforminv_nsym(S_rs,S_HSSinv,S_FTinv);

elseif strcmp(opts,'d2HSS')
    
    X1 = X(I1_sk,:); X2 = X(I2_sk,:);
    m1 = size(X1,1); 
    
    % If opts == 'd2HSS', then FI is built using the Schur complement
    % formulae. It is assumed that E1 and E2 are dense matrices. 
    % The diagonal blocks of F are computed as E_dg + K_offd. 
    
    acc = params.acc;
    JE1 = params.JE1; JE2 = params.JE2; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % (I) Computing F_sk and F_rs in HSS form: 
    
    % (0) SORT: we sort points and separate skeleton and residual sets.
    Bl = false(size(J)); Bl(J(1:k)) = true; 
    
    B1 = Bl(1:m1); B2 = Bl(m1+1:end); 
    J1 = J(J<=m1); J2 = J(J>m1) - m1; 
    B1 = B1(J1); B2 = B2(J2); 
    X1_s = X1(J1,:); X2_s = X2(J2,:); 
    
    E1 = E1(J1,J1); E2 = E2(J2,J2); 
    
    % Consolidate points
    Xrs  = [X1_s(~B1,:) ; X2_s(~B2,:)];
    
    [~,Jrs] = sort_residual(Xrs,lev+1);
    
    J1 = JE1(J1); J2 = JE2(J2); 
    
    JFr = [J1(~B1) (J2(~B2) + length(B1))];
    FI.J = [J1(B1) (J2(B2) + length(B1)) JFr(Jrs)]; 
    
    Tsk = LOCAL_rs_tree(sum(B1),sum(B2));
    Trs = LOCAL_rs_tree(sum(~B1),sum(~B2));
    
    % (3) OFFD
    %F_sk  = HSS1D_compress_fsym_offd_giventree(Xsk,Tsk,acc,'lap',params,1); 
    %F_sk{BOX.M_SELF,2} = E1(B1,B1);   F_sk{BOX.M_SELF,3} = E2(B2,B2); 
    Fd_sk = [E1(B1,B1) Kernel_Eval(X1_s(B1,:),X2_s(B2,:),params) ; Kernel_Eval(X2_s(B2,:),X1_s(B1,:),params) E2(B2,B2)];
    
    Fd_rs = [E1(~B1,~B1) Kernel_Eval(X1_s(~B1,:),X2_s(~B2,:),params) ;...
        Kernel_Eval(X2_s(~B2,:),X1_s(~B1,:),params) E2(~B2,~B2)];
    Fd_rs = Fd_rs(Jrs,Jrs);
    
    %F_rs  = HSS1D_compress_brute_nsym_giventree(Fd_rs,Trs,acc);
    F_rs   = HSS1D_compress_nsym('brute',Fd_rs,Trs,acc); 
    %F_sk  = HSS1D_compress_brute_nsym_giventree(Fd_sk,Tsk,acc);
    F_sk   = HSS1D_compress_nsym('brute',Fd_sk,Tsk,acc); 
    
    % Invert F_rs using HSS1D inversion
    [Frs_HSSinv,Frs_FTinv] = HSS1D_invert_nsym(F_rs,acc);
    FI.rs_Inv = HSS1D_transforminv_nsym(F_rs,Frs_HSSinv,Frs_FTinv);    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (II) Off-diagonal blocks of F as low rank using ID_rand
    
    Fr2s  = [E1(B1,~B1) Kernel_Eval(X1_s(B1,:),X2_s(~B2,:),params) ; Kernel_Eval(X2_s(B2,:),X1_s(~B1,:),params) E2(B2,~B2)];
    Fs2r  = [E1(~B1,B1) Kernel_Eval(X1_s(~B1,:),X2_s(B2,:),params) ; Kernel_Eval(X2_s(~B2,:),X1_s(B1,:),params) E2(~B2,B2)];
    
    T_r2s = LRID(Fr2s,acc);
    T_s2r = LRID(Fs2r,acc);
    
    T_r2s.V = T_r2s.V(Jrs,:); T_s2r.U = T_s2r.U(Jrs,:); 
    
    FI.rs2sk = T_r2s; 
    FI.sk2rs = T_s2r;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (III) Schur Complement Matrix S_rs
    
    U2 = HSS1D_apply_nsym(FI.rs_Inv,T_s2r.U); 
    U = -LR_apply(T_r2s,U2); 
    V = T_s2r.V; 
    
    M = HSS1D_compressLR_nsym(U,V,F_sk,acc);
    S_rs = HSS1D_sum_nsym(F_sk,M,acc); 
    S_rs = HSS1D_recompress_nsym(S_rs,acc);
 
    [S_HSSinv,S_FTinv] = HSS1D_invert_nsym(S_rs,acc);
    FI.SrsInv = HSS1D_transforminv_nsym(S_rs,S_HSSinv,S_FTinv);

elseif strcmp(opts,'HSStest')
    X1 = X(I1_sk,:); X2 = X(I2_sk,:);
    
    % If opts == 'HSS', then FI is built using the Schur complement
    % formulae. It is assumed that E1 and E2 are HSS matrices. 
    % The diagonal blocks of F are computed as E_dg + K_offd using HSS1DSUM. 
    
    acc = params.acc; 
    JE1 = params.JE1; JE2 = params.JE2; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % (I) Computing F_sk and F_rs in HSS form: 
    X1_s = X1(JE1,:); X2_s = X2(JE2,:); 
    
    % (0) SORT: we sort points and separate skeleton and residual sets.
    Bl = false(size(J')); Bl(J(1:k)) = true; 
    
    m1 = size(X1,1); 
    B1 = Bl(1:m1); B2 = Bl(m1+1:end); 
    B1 = B1(JE1); B2 = B2(JE2); 
        
    E1d = HSS1D_apply_nsym(E1,eye(E1{BOX.END2,1}));
    E2d = HSS1D_apply_nsym(E2,eye(E2{BOX.END2,1}));
        
    % (1) SPLIT: we split the E HSS matrices into skeleton and residual sets
    [E1_sk,E1_rs] = HSS1D_split_nsym(E1,B1,B1,acc); 
    [E2_sk,E2_rs] = HSS1D_split_nsym(E2,B2,B2,acc); 
    
    k1 = sum(B1); k2 = sum(B2); k = k1+k2; 
    r1 = sum(~B1); r2 = sum(~B2); r = r1+r2; 
    
    V_sk = randn(k,10); V_sk = V_sk/norm(V_sk); 
    V_rs = randn(r,10); V_rs = V_rs/norm(V_rs); 
    
    fprintf('\nSplit test\n')
    norm(E1d(B1,B1)*V_sk(1:k1,:) - HSS1D_apply_nsym(E1_sk,V_sk(1:k1,:)))
    norm(E1d(~B1,~B1)*V_rs(1:r1,:) - HSS1D_apply_nsym(E1_rs,V_rs(1:r1,:)))
    
    % Consolidate points
    Xsk  = [X1_s(B1,:) ; X2_s(B2,:)]; 
    Xrs  = [X1_s(~B1,:) ; X2_s(~B2,:)];  
    
    % (2) MERGE
    E_sk              = HSS1D_merge_nsym(E1_sk,E2_sk,Xsk,lev,acc,'concatenate');
    [E_rs,Xrs_s,Jrs]  = HSS1D_merge_nsym(E1_rs,E2_rs,Xrs,lev,acc,'align');
    
    fprintf('\nMerge test\n')
    norm([E1d(B1,B1)*V_sk(1:k1,:) ; E2d(B2,B2)*V_sk(k1+1:end,:)] - HSS1D_apply_nsym(E_sk,V_sk)) 
    Esk_dg = [E1d(B1,B1) zeros(k1,k2) ; zeros(k2,k1) E2d(B2,B2)]; 
    Ers_dg = [E1d(~B1,~B1) zeros(r1,r2) ; zeros(r2,r1) E2d(~B2,~B2)]; 
    Yrs = Ers_dg*V_rs;
    YHrs(Jrs,:) = HSS1D_apply_nsym(E_rs,V_rs(Jrs,:));
    norm(Yrs - YHrs)
    
    Brs = [true(1,sum(~B1)) false(1,sum(~B2))]; 
    Brs = Brs(Jrs); 
    
    % (3) OFFD and SUM
    %Ksk_off  = HSS1D_compress_nsym_offd_giventree(Xsk,E_sk,acc,'lap',params,1); 
    Ksk_off   = HSS1D_compress_nsym('green',Xsk,E_sk,acc,params,'offd'); 
    Krs_off  = HSS1D_compress_nsym_offd_leaf(Xrs_s,E_rs,Brs,acc,params);
    
    fprintf('\nOffd test\n')
    Kskod = [zeros(k1,k1) Kernel_Eval(Xsk(1:k1,:),Xsk(k1+1:end,:),params) ; Kernel_Eval(Xsk(k1+1:end,:),Xsk(1:k1,:),params) zeros(k2,k2)];
    Krsod = [zeros(r1,r1) Kernel_Eval(Xrs(1:r1,:),Xrs(r1+1:end,:),params) ; Kernel_Eval(Xrs(r1+1:end,:),Xrs(1:r1,:),params) zeros(r2,r2)];
    
    norm(Kskod*V_sk - HSS1D_apply_nsym(Ksk_off,V_sk)) 
    norm(Krsod(Jrs,Jrs)*V_rs - HSS1D_apply_nsym(Krs_off,V_rs))
       
    F_sk  = HSS1D_sum_nsym(E_sk,Ksk_off,acc); 
    F_rs  = HSS1D_sum_nsym(E_rs,Krs_off,acc);
    F_rs = HSS1D_recompress_nsym(F_rs,acc);  
    
    fprintf('\nSum test\n')
    Fd_sk = (Kskod + Esk_dg);
    Fd_rs = (Krsod(Jrs,Jrs) + Ers_dg(Jrs,Jrs));
    norm(Fd_sk*V_sk - HSS1D_apply_nsym(F_sk,V_sk))
    norm(Fd_rs*V_rs - HSS1D_apply_nsym(F_rs,V_rs))
   
    [Frs_HSSinv,Frs_FTinv] = HSS1D_invert_nsym(F_rs,acc);
    FI.rs_Inv = HSS1D_transforminv_nsym(F_rs,Frs_HSSinv,Frs_FTinv);    
    
    fprintf('\nFrs_Inv test\n')
    norm(Fd_rs\V_rs - HSS1D_applyinv_nsym(F_rs,Frs_HSSinv,Frs_FTinv,V_rs))/norm(Fd_rs\V_rs)
    norm(Fd_rs\V_rs - HSS1D_apply_nsym(FI.rs_Inv,V_rs))/norm(Fd_rs\V_rs)
    
    J1 = JE1; J2 = JE2; 
    JFr = [J1(~B1) (J2(~B2) + length(B1))];
    FI.J = [J1(B1) (J2(B2) + length(B1)) JFr(Jrs)];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (II) Off-diagonal blocks of F as low rank using ID_rand
    
    Fr2s  = [E1d(B1,~B1) Kernel_Eval(X1_s(B1,:),X2_s(~B2,:),params) ; Kernel_Eval(X2_s(B2,:),X1_s(~B1,:),params) E2d(B2,~B2)];
    Fs2r  = [E1d(~B1,B1) Kernel_Eval(X1_s(~B1,:),X2_s(B2,:),params) ; Kernel_Eval(X2_s(~B2,:),X1_s(B1,:),params) E2d(~B2,B2)];
    
    Td_r2s = LRID(Fr2s,acc);
    Td_s2r = LRID(Fs2r,acc);
    
     Td_r2s.V = Td_r2s.V(Jrs,:); Td_s2r.U = Td_s2r.U(Jrs,:); 
    
    [T_r2s,T_s2r] = LOCAL_LR_Tr2s(E1,E2,B1,B2,X1_s,X2_s,lev,acc,params);
    
    T_r2s.V = T_r2s.V(Jrs,:); T_s2r.U = T_s2r.U(Jrs,:); 
    
    FI.rs2sk = T_r2s;
    FI.sk2rs = T_s2r; 
    
    fprintf('\n Error Fr2s / Fs2r LR \n')
    norm(T_r2s.U*T_r2s.V.' - Fr2s(:,Jrs))/norm(Fr2s)
    norm(T_s2r.U*T_s2r.V.' - Fs2r(Jrs,:))/norm(Fs2r)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (III) Schur Complement Matrix S_rs
    
    U2 = HSS1D_apply_nsym(FI.rs_Inv,T_s2r.U); 
    U = -LR_apply(T_r2s,U2); 
    V = T_s2r.V; 
    
    %S_rs = HSS1D_LRUpdate_nsym(F_sk,U,V,acc); 
    M = HSS1D_compressLR_nsym(U,V,F_sk,acc);
    S_rs = HSS1D_sum_nsym(F_sk,M,acc); 
    S_rs = HSS1D_recompress_nsym(S_rs,acc); 
    
    [S_HSSinv,S_FTinv] = HSS1D_invert_nsym(S_rs,acc);
    FI.SrsInv = HSS1D_transforminv_nsym(S_rs,S_HSSinv,S_FTinv);
    
    fprintf('S_rs test')
    Sd_rs = Fd_sk + U*V.'; 
    norm(Sd_rs*V_sk - HSS1D_apply_nsym(S_rs,V_sk))/norm(Sd_rs*V_sk)
    norm(U*(V.'*V_sk) - HSS1D_apply_nsym(M,V_sk))/norm(U*(V.'*V_sk))
    
    fprintf('\nSrsInv test\n')
    norm(Sd_rs\V_sk - HSS1D_applyinv_nsym(S_rs,S_HSSinv,S_FTinv,V_sk))/norm(Sd_rs\V_sk)
    norm(Sd_rs\V_sk - HSS1D_apply_nsym(FI.SrsInv,V_sk))/norm(Sd_rs\V_sk)
else
    X1 = X(I1_sk,:); X2 = X(I2_sk,:);
    
    % If opts == 'blocks', then FI is built densely using the Schur complement
    % formulae. It is assumed that E1 and E2 are dense matrices. 
    
    acc = params.acc; 
    JE1 = params.JE1; JE2 = params.JE2; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % (I) Computing F_sk and F_rs: 
     X1_s = X1(JE1,:); X2_s = X2(JE2,:); 
     m1 = size(X1,1); 
    
    if strcmp(opts,'block')
        % (0) SORT: we sort points and separate skeleton and residual sets
        Bl = false(size(J')); Bl(J(1:k)) = true; 

        B1 = Bl(1:m1); B2 = Bl(m1+1:end); 
        B1 = B1(JE1); B2 = B2(JE2); 
   
        J1 = JE1; J2 = JE2;
    else
        Bl = false(size(J)); Bl(J(1:k)) = true; 
    
        B1 = Bl(1:m1); B2 = Bl(m1+1:end); 
        J1 = J(J<=m1); J2 = J(J>m1) - m1; 
        B1 = B1(J1); B2 = B2(J2); 
        X1_s = X1(J1,:); X2_s = X2(J2,:); 
    
        E1 = E1(J1,J1); E2 = E2(J2,J2); 
    
        J1 = JE1(J1); J2 = JE2(J2);     
    end
    
    Xrs  = [X1_s(~B1,:) ; X2_s(~B2,:)];
    [~,Jrs] = sort_residual(Xrs,lev+1);
    JFr = [J1(~B1) (J2(~B2) + length(B1))];
    FI.J = [J1(B1) (J2(B2) + length(B1)) JFr(Jrs)];  
    
    F_sk  = [E1(B1,B1) Kernel_Eval(X1_s(B1,:),X2_s(B2,:),params) ; Kernel_Eval(X2_s(B2,:),X1_s(B1,:),params) E2(B2,B2)];
    F_rs  = [E1(~B1,~B1) Kernel_Eval(X1_s(~B1,:),X2_s(~B2,:),params) ; Kernel_Eval(X2_s(~B2,:),X1_s(~B1,:),params) E2(~B2,~B2)];
    
    % Invert F_rs
    F_rs = F_rs(Jrs,Jrs); 
    FI.rs_Inv = LOCAL_Invert_sym(F_rs); 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (II) Off-diagonal blocks of F as low rank using ID_rand
    
    Fr2s  = [E1(B1,~B1) Kernel_Eval(X1_s(B1,:),X2_s(~B2,:),params) ; Kernel_Eval(X2_s(B2,:),X1_s(~B1,:),params) E2(B2,~B2)];
    Fs2r  = [E1(~B1,B1) Kernel_Eval(X1_s(~B1,:),X2_s(B2,:),params) ; Kernel_Eval(X2_s(~B2,:),X1_s(B1,:),params) E2(~B2,B2)];
    
    T_r2s = LRID(Fr2s,acc);
    T_s2r = LRID(Fs2r,acc);
    
     T_r2s.V = T_r2s.V(Jrs,:); T_s2r.U = T_s2r.U(Jrs,:); 
    
    FI.rs2sk = T_r2s; 
    FI.sk2rs = T_s2r;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (III) Schur Complement Matrix S_rs
    
    % (III) Schur Complement matrix
    U2 = FI.rs_Inv*T_s2r.U; 
    U = -LR_apply(T_r2s,U2); 
    V = T_s2r.V;
    
    S_rs = F_sk + U*V.'; 
    
    % Inverse
    FI.SrsInv = inv(S_rs);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TREE = LOCAL_sk_tree(T1,T2)
global BOX 

n1 = size(T1,2);
n2 = size(T2,2); 
ntot = T1{BOX.END2,1} + T2{BOX.END2,1}; 
nboxes = n1+n2+1;
TREE = cell(50,nboxes);

% Construct the top node.
TREE{BOX.LEVEL,1} = 0;
TREE{BOX.PARENT,1} = NaN;
TREE{BOX.C1,1} = 2;
TREE{BOX.C2,1} = 3;
TREE{BOX.END1,1} = 1;
TREE{BOX.END2,1} = ntot; 

% roots of T1 and T2 are now children of new root node
T1{BOX.PARENT,1} = 0; 
T2{BOX.PARENT,1} = 0; 

numlev1 = zeros(1,T1{BOX.LEVEL,n1}+1); 
numlev2 = zeros(1,T2{BOX.LEVEL,n2}+1); 

depth = max(T1{BOX.LEVEL,n1},T2{BOX.LEVEL,n2}); 

for i=1:n1
   numlev1(T1{BOX.LEVEL,i}+1) = numlev1(T1{BOX.LEVEL,i}+1) + 1;  
   T1{BOX.LEVEL,i} = T1{BOX.LEVEL,i} + 1; 
end

for i=1:n2
   numlev2(T2{BOX.LEVEL,i}+1) = numlev2(T2{BOX.LEVEL,i}+1 ) + 1;
   T2{BOX.LEVEL,i} = T2{BOX.LEVEL,i} + 1; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Concatenate TREE = [T1 0 ; 0 T2] 

    n_j1 = 0; n_j2 = 0; n_jtot = 1; 
    for li = 1:depth+1
       % Modify info for T1 boxes
       i1_loc = (n_j1+1:n_j1+numlev1(li));
       
       for j = 1:numlev1(li)
          % Parent
          T1{BOX.PARENT,i1_loc(j)} = T1{BOX.PARENT,i1_loc(j)} + n_jtot; 
          
          if li>1
              T1{BOX.PARENT,i1_loc(j)} = T1{BOX.PARENT,i1_loc(j)} - n_j2 - numlev1(li-1); 
          end
          
          % Children (if the box is not a leaf)
          if li<T1{BOX.LEVEL,n1}
              if (T1{BOX.C1,i1_loc(j)}>0)
                 T1{BOX.C1,i1_loc(j)} = T1{BOX.C1,i1_loc(j)} + n_j2 + 1 + numlev2(li);
              end
              if (T1{BOX.C2,i1_loc(j)}>0)
                 T1{BOX.C2,i1_loc(j)} = T1{BOX.C2,i1_loc(j)} + n_j2 + 1 + numlev2(li); 
              end
          end
       end
       
       % Add T1 boxes to E on level li
       TREE(:,n_jtot+1:n_jtot+numlev1(li)) = T1(:,i1_loc); 
       n_j1   = n_j1  + numlev1(li); 
       n_jtot = n_jtot+ numlev1(li);
       
       % Modify info for T2 boxes
       i2_loc = (n_j2+1:n_j2+numlev2(li)); 
       for j = 1:numlev2(li)
           % Parent
           T2{BOX.PARENT,i2_loc(j)} = T2{BOX.PARENT,i2_loc(j)} + n_jtot; 
           
           T2{BOX.PARENT,i2_loc(j)} = T2{BOX.PARENT,i2_loc(j)} - n_j1; 
          
           % Children (if the box is not a leaf)
           if li<T2{BOX.LEVEL,n2}
               if (T2{BOX.C1,i2_loc(j)}>0)
                  T2{BOX.C1,i2_loc(j)} = T2{BOX.C1,i2_loc(j)} + n_j1 + 1 + numlev1(li+1); 
               end
               if (T2{BOX.C2,i2_loc(j)}>0)
                  T2{BOX.C2,i2_loc(j)} = T2{BOX.C2,i2_loc(j)} + n_j1 + 1 + numlev1(li+1); 
               end
           end
           T2{BOX.END1,i2_loc(j)}  = T2{BOX.END1,i2_loc(j)}  + T1{BOX.END2,1}; % Update starting index
           T2{BOX.I_SKUP,i2_loc(j)} = T2{BOX.I_SKUP,i2_loc(j)} + T1{BOX.END2,1}; % indskel_out
           T2{BOX.I_SKDN,i2_loc(j)} = T2{BOX.I_SKDN,i2_loc(j)} + T1{BOX.END2,1}; % indskel_in
       end
       
       % Add T2 boxes to E on level li
       TREE(:,n_jtot+1:n_jtot+numlev2(li)) = T2(:,i2_loc); 
       n_j2   = n_j2  + numlev2(li); 
       n_jtot = n_jtot+ numlev2(li);
       
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TREE = LOCAL_rs_tree(n1,n2)
global BOX 

ntot = n1+n2; 
TREE = cell(50,3);

% Construct the top node.
TREE{BOX.LEVEL,1} = 0;
TREE{BOX.PARENT,1} = NaN;
TREE{BOX.C1,1} = 2;
TREE{BOX.C2,1} = 3;
TREE{BOX.END1,1} = 1;
TREE{BOX.END2,1} = ntot; 

% Nodes 2 and 3
TREE{BOX.LEVEL,2} = 1; TREE{BOX.LEVEL,3} = 1; 
TREE{BOX.PARENT,2} = 1; TREE{BOX.PARENT,3} = 1;  
TREE{BOX.END1,2} = 1; TREE{BOX.END2,2} = n1; 
TREE{BOX.END1,3} = n1+1; TREE{BOX.END2,3} = n2; 

TREE{BOX.C1,2} = -1; TREE{BOX.C2,2} = -1;
TREE{BOX.C1,3} = -1; TREE{BOX.C2,3} = -1;

% Split in half
%nbox_max = ceil(min(n1,n2)/2);
nbox_max = 56; 
% parameters to keep subdividing the tree, if necessary
ibox_last = 1;
ncreated  = 2;
ilevel    = 1;

% Create smaller nodes via hierarchical subdivision.
% We sweep one level at a time.
while (ncreated > 0)
  ibox_first = ibox_last + 1;
  ibox_last  = ibox_last + ncreated;
  ncreated   = 0;
  for ibox = ibox_first:ibox_last
    nbox = TREE{BOX.END2,ibox};
    if (nbox > nbox_max)
      nhalf             = ceil(nbox/2);
      ibox_son1         = ibox_last + ncreated + 1;
      ibox_son2         = ibox_last + ncreated + 2;
      TREE{BOX.C1,ibox}      = ibox_son1;
      TREE{BOX.C2,ibox}      = ibox_son2;
      TREE{BOX.LEVEL,ibox_son1} = ilevel+1;
      TREE{BOX.LEVEL,ibox_son2} = ilevel+1;
      TREE{BOX.PARENT,ibox_son1} = ibox;
      TREE{BOX.PARENT,ibox_son2} = ibox;
      TREE{BOX.C1,ibox_son1} = -1;
      TREE{BOX.C1,ibox_son2} = -1;
      TREE{BOX.C2,ibox_son1} = -1;
      TREE{BOX.C2,ibox_son2} = -1;
      TREE{BOX.END1,ibox_son1} = TREE{BOX.END1,ibox};
      TREE{BOX.END1,ibox_son2} = TREE{BOX.END1,ibox} + nhalf;
      TREE{BOX.END2,ibox_son1} = nhalf;
      TREE{BOX.END2,ibox_son2} = TREE{BOX.END2,ibox} - nhalf;
      ncreated          = ncreated + 2;
    end
  end
  ilevel = ilevel + 1;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [T_r2s,T_s2r] = LOCAL_LR_Tr2s(E1,E2,B1,B2,X1_s,X2_s,lev,acc,params)

Br2s_dn = cell(2,2); Br2s_up = Br2s_dn; 
Bs2r_dn = Br2s_dn; Bs2r_up = Bs2r_dn; 

% Bool arrays for restricted tree traversal 
[Br2s_dn{1,1},Br2s_up{1,1}] = HSS1D_apply_submatrix_boxes(E1,B1,~B1);
[Br2s_dn{2,2},Br2s_up{2,2}] = HSS1D_apply_submatrix_boxes(E2,B2,~B2);  
[Bs2r_dn{1,1},Bs2r_up{1,1}] = HSS1D_apply_submatrix_boxes(E1,~B1,B1);
[Bs2r_dn{2,2},Bs2r_up{2,2}] = HSS1D_apply_submatrix_boxes(E2,~B2,B2);  

k1 = sum(B1); k2 = sum(B2); r1 = sum(~B1); r2 = sum(~B2); 
Bl{1,1} = B1; Bl{2,2} = B2; 

% Restricted Compress of Off-diagonal blocks
if (k1 == k2 && r1 == r2)
    Bl{1,2} = B1; Bl{2,1} = B2; 
    
    % Joint vectors X12 = [I1_s,I2_r] and X21 = [I2_s,I1_r]
    X12(B1,:)=X1_s(B1,:); X12(~B1,:)=X2_s(~B2,:);
    X21(B2,:)=X2_s(B2,:); X21(~B2,:)=X1_s(~B1,:);
    % Bool arrays for restricted tree traversal 
    Br2s_dn{1,2} = Br2s_dn{1,1}; Br2s_dn{2,1} = Br2s_dn{2,2}; 
    Br2s_up{1,2} = Br2s_up{1,1}; Br2s_up{2,1} = Br2s_up{2,2}; 
    Bs2r_dn{1,2} = Bs2r_dn{2,2}; Bs2r_dn{2,1} = Bs2r_dn{1,1}; 
    Bs2r_up{1,2} = Bs2r_up{2,2}; Bs2r_up{2,1} = Bs2r_up{1,1}; 

    % K[I1_s,I2_r] and K[I2_s,I1_r]
    K12_r2s = HSS1D_compress_submatrix_nsym_giventree(X12,E1,Br2s_dn{1,1},Br2s_up{1,1},acc/100,params);
    K21_r2s = HSS1D_compress_submatrix_nsym_giventree(X21,E2,Br2s_dn{2,2},Br2s_up{2,2},acc/100,params);
    
    % K[I1_r,I2_s] and K[I2_r,I1_s]
    K12_s2r = HSS1D_compress_submatrix_nsym_giventree(X21,E2,Bs2r_dn{1,2},Bs2r_up{1,2},acc/100,params);
    K21_s2r = HSS1D_compress_submatrix_nsym_giventree(X12,E1,Bs2r_dn{2,1},Bs2r_up{2,1},acc/100,params);     
else
    
    Bl{1,2} = [true(1,k1) false(1,r2)]; 
    Bl{2,1} = [true(1,k2) false(1,r1)];
    % Joint vectors X12 = [I1_s,I2_r] and X21 = [I2_s,I1_r]
    X12 = [X1_s(B1,:) ; X2_s(~B2,:)]; T12 = LOCAL_rs_tree(k1,r2); 
    X21 = [X2_s(B2,:) ; X1_s(~B1,:)]; T21 = LOCAL_rs_tree(k2,r1); 
    % Bool arrays for restricted tree traversal 
    [Br2s_dn{1,2},Br2s_up{1,2}] = HSS1D_apply_submatrix_boxes(T12,Bl{1,2},~Bl{1,2});
    [Br2s_dn{2,1},Br2s_up{2,1}] = HSS1D_apply_submatrix_boxes(T21,Bl{2,1},~Bl{2,1});
    
    [Bs2r_dn{1,2},Bs2r_up{1,2}] = HSS1D_apply_submatrix_boxes(T21,~Bl{2,1},Bl{2,1});
    [Bs2r_dn{2,1},Bs2r_up{2,1}] = HSS1D_apply_submatrix_boxes(T12,~Bl{1,2},Bl{1,2});
    
    % K[I1_s,I2_r] and K[I2_s,I1_r]
    K12_r2s = HSS1D_compress_submatrix_nsym_giventree(X12,T12,Br2s_dn{1,2},Br2s_up{1,2},acc/100,params);    
    K21_r2s = HSS1D_compress_submatrix_nsym_giventree(X21,T21,Br2s_dn{2,1},Br2s_up{2,1},acc/100,params);
    
    % K[I1_r,I2_s] and K[I2_r,I1_s]
    K12_s2r = HSS1D_compress_submatrix_nsym_giventree(X21,T21,Bs2r_dn{1,2},Bs2r_up{1,2},acc/100,params);
    K21_s2r = HSS1D_compress_submatrix_nsym_giventree(X12,T12,Bs2r_dn{2,1},Bs2r_up{2,1},acc/100,params);    
end

Blt = cell(size(Bl)); Bl_up = Blt; Blt_up = Blt; 
for i=1:2
    for j=1:2
        Blt{i,j}    = ~Bl{j,i};
        Bl_up{i,j}  = ~Bl{i,j};
        Blt_up{i,j} = Bl{j,i};
    end
end

% Transpose matrices
E1t = HSS1D_transpose(E1); 
E2t = HSS1D_transpose(E2); 
K12t_r2s = HSS1D_transpose(K21_r2s); 
K21t_r2s = HSS1D_transpose(K12_r2s); 
K12t_s2r = HSS1D_transpose(K21_s2r); 
K21t_s2r = HSS1D_transpose(K12_s2r); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Block HSS (BHSS) function handle applies
mvBr2s   = @(x) BHSS1D_apply_nsym(E1,K12_r2s,K21_r2s,E2,Br2s_dn,Br2s_up,Bl,Bl_up,x);

tmp = Br2s_up{1,2}; Br2s_up{1,2} = Br2s_up{2,1}; Br2s_up{2,1} = tmp; 
tmp = Br2s_dn{1,2}; Br2s_dn{1,2} = Br2s_dn{2,1}; Br2s_dn{2,1} = tmp; 

mvBr2st = @(x) BHSS1D_apply_nsym(E1t,K12t_r2s,K21t_r2s,E2t,Br2s_up,Br2s_dn,Blt,Blt_up,x);

mvBs2r = @(x) BHSS1D_apply_nsym(E1,K12_s2r,K21_s2r,E2,Bs2r_dn,Bs2r_up,Blt,Blt_up,x);

tmp = Bs2r_up{1,2}; Bs2r_up{1,2} = Bs2r_up{2,1}; Bs2r_up{2,1} = tmp; 
tmp = Bs2r_dn{1,2}; Bs2r_dn{1,2} = Bs2r_dn{2,1}; Bs2r_dn{2,1} = tmp; 

mvBs2rt = @(x) BHSS1D_apply_nsym(E1t,K12t_s2r,K21t_s2r,E2t,Bs2r_up,Bs2r_dn,Bl,Bl_up,x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Randomized ID parameters
    
ID_params.form = 'BHSS'; 
ID_params.acc=acc/10;
dist = 'norm'; 
m1 = k1 + k2;
n1 = r1 + r2;
rho = params.rho; 

% Rank guess q and stride C
if rho(lev+1)>0
    ID_params.q = rho(lev+1); 
    ID_params.C = 50; 
elseif sum(rho>0)<3
    [ID_params.q,ID_params.C] = predict_rank(n1,acc);
else
    ID_params.q = sqrt(2)*rho(lev+3);
    ID_params.C = ceil((rho(lev+2) - rho(lev+3))/2) + 10;
end

% (par == acc) -> Rand ID finds rank up to that accuracy
% (par == q  ) -> Rand ID for fixed rank q
ID_params.par = acc; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rand ID to obtain T_r2s and T_s2r in Low Rank form

T_r2s = LRID_rand(m1,n1,mvBr2s,mvBr2st,dist,ID_params);

ID_params.q = max(ID_params.q,size(T_r2s.U,2)); 

T_s2r = LRID_rand(n1,m1,mvBs2r,mvBs2rt,dist,ID_params);

end