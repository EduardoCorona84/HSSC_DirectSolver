function [E,Js_J2F] = HSS2D_build_E_fsym(FI,J,T,k,acc,opts)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%     FUNCTION CALL:
%         [E,Js_J2F] = HSS2D_build_E_fsym(FI,J,T,k,acc,opts)
%
%     DESCRIPTION:
%         This function is an implementation of Algorithm 5 in the paper - which corresponds to
%         line 10 of Algorithm 2. For a given box, it computes E as a function of the
%         interpolation operators and the blocks of F^-1. It assumes that the associated kernel
%         is symmetric (fsym).
%
%     INPUT:
%         FI      <nxn or struct>     (output of HSS2D_build_Finv_fsym.m) The blocks of F^-1
%                                         needed to compute E.
%         J       <kx1 int>           [I_sk I_rs] Permutation vector that places skeleton points
%                                         first, residual points second.
%         T       <kxk or lr>         Interpolation matrix (stored either densely or as a
%                                         low-rank (U,V) pair).
%         k       <int>               Size of skeleton.
%         acc     <float>             Desired accuracy.
%
%         opts    <struct>            .mat determines how E is stored. If 'dense' or 'block', E
%                                         is dense and .mat.LR specifies if the interpolation
%                                         operations are stored as lowrank. If .mat = 'HSS', E
%                                         is HSS1D and interpolation operators are assumed to be
%                                         lowrank.
%
%     OUTPUT:
%         E   <kxk or HSS1D>   Schur Complement matrix used to represent A^-1.
%
%                         EInv = R*FInv*L = [I T] * FI * [I ; T.']
%                         E = inv(EInv)
%




if strcmp(opts.mat,'dense')
    
    if (opts.LR == 0)
        R(:,J) = [ eye(k) T.U]; 
        L = R.'; 
        EI = R*(FI*L);
        
        % Direct inversion
        E = inv(EI); 
    else
        Phi = FI(J,J);
        U_cell = cell(3,1); V_cell = cell(3,1); 
        
        V_cell{1} = Phi(k+1:end,1:k).'*T.V;  U_cell{1} = T.U;
        V_cell{2} = T.U;                    U_cell{2} = Phi(1:k,k+1:end)*(T.V); 
        V_cell{3} = T.U;                    
        
        U_cell{3} = T.U*(T.V.')*Phi(k+1:end,k+1:end)*T.V; 
        
        % Inversion 
        EI = Phi(1:k,1:k) + U_cell{1}*V_cell{1}.' + U_cell{2}*V_cell{2}.' + U_cell{3}*V_cell{3}.'; 
        E = inv(EI); 
        
        %{
        SMW (Sherman-Morrison-Woodbury)
        S = inv(Phi(1:k,1:k)); 
        m = size(U,2); 
        MI = (eye(m) + V'*(Phi(1:k,1:k)\U)); 
        E = S*(eye(k) - U*(MI\V')*S); 
        %}
    end
    
    Js_J2F = 1:k; 
    
elseif strcmp(opts.mat,'HSS')
    
        % phi_rs = FI.rs_Inv*(eye(length(FI.rs_Inv)) + FI.sk2rs*FI.SrsInv*FI.rs2sk*FI.rs_Inv);
        % phi_rs2sk = -FI.SrsInv*FI.rs2sk*FI.rs_Inv; 
        % phi_sk2rs = -FI.rs_Inv*FI.sk2rs*FI.SrsInv; 
        
        %fprintf('\n E matrix as LR Update of SrsInv \n')
        J_F = FI.J; 
        Ind = 1:length(J); Jinv(J) = Ind; 
        
        J_J2F = Jinv(J_F);
        Js_J2F = J_J2F(1:k); 
        Jr_J2F = J_J2F(k+1:end) - k; 
        
        % Permute T to make indices the same as in F
        T.U = T.U(Js_J2F,:); 
        T.V = T.V(Jr_J2F,:); 
        
        U_cell = cell(3,1); V_cell = cell(3,1); 
        
        Vtmp1 = HSS1D_apply_fsym(FI.rs_Inv,T.V); 
        Vtmp1 = LR_apply(FI.rs2sk,Vtmp1); 
        
        V_cell{1} = HSS1D_apply_fsym(FI.SrsInv,Vtmp1); 
        U_cell{1} = -T.U;
        
        V_cell{2} = U_cell{1};          
        U_cell{2} = V_cell{1}; 
        
        V_cell{3} = T.U;  
        
        Fs2r.U = FI.rs2sk.V; Fs2r.V = FI.rs2sk.U; 
        Utmp3 = T.V + LR_apply(Fs2r,U_cell{2});
        Utmp3 = HSS1D_apply_fsym(FI.rs_Inv,Utmp3); 
        U_cell{3} = LR_apply(T,Utmp3); 
        
        [U,V] = LR_recomp(U_cell,V_cell,acc/10);
        M = HSS1D_compressLR_fsym(U,V,FI.SrsInv,acc);
        M = HSS1D_symmetrize_fsym(M); 
        
        EI = HSS1D_sum_fsym(FI.SrsInv,M,acc); 
        EI = HSS1D_recompress_fsym(EI,acc);
        EI = HSS1D_fsym_to_nsym(EI); 
        
        [E_HSSinv,E_FTinv] = HSS1D_invert_nsym(EI,acc);
        E = HSS1D_transforminv_nsym(EI,E_HSSinv,E_FTinv);       
        E = HSS1D_nsym_to_fsym(E,acc); 
    
elseif strcmp(opts.mat,'HSStest')
    
    fprintf('\n HSStest: BuildE \n')
    Vsk = randn(k,10); Vsk = Vsk./norm(Vsk); 
    
    %fprintf('\n E matrix as LR Update of SrsInv \n')
    J_F = FI.J; 
    Ind = 1:length(J); Jinv(J) = Ind; 
        
    J_J2F = Jinv(J_F);
    Js_J2F = J_J2F(1:k); 
    Jr_J2F = J_J2F(k+1:end) - k; 
        
    % Permute T to make indices the same as in F
    T.U = T.U(Js_J2F,:); 
    T.V = T.V(Jr_J2F,:); 
    
    %Convert HSS matrices to dense for comparison
    FId_rs_Inv = HSS1D_dense_matrix(FI.rs_Inv,1); 
    FId_SrsInv = HSS1D_dense_matrix(FI.SrsInv,1); 
        
    U_cell = cell(3,1); V_cell = cell(3,1); 
    
    % U1*V1.'
    Vtmp1 = HSS1D_apply_fsym(FI.rs_Inv,T.V); 
    Vtmp1 = LR_apply(FI.rs2sk,Vtmp1); 
        
    V_cell{1} = HSS1D_apply_fsym(FI.SrsInv,Vtmp1); 
    U_cell{1} = -T.U;
    
    V1d = FId_SrsInv*(LR_apply(FI.rs2sk,FId_rs_Inv*T.V)); 
    U1d = -T.U; 
    
    % U2 = V1 and V2 = U1
    V_cell{2} = U_cell{1};          
    U_cell{2} = V_cell{1}; 
    
    % U3*V3.'
    V_cell{3} = T.U;  
    V3d = T.U; 
        
    Fs2r.U = FI.rs2sk.V; Fs2r.V = FI.rs2sk.U; 
    Utmp3 = T.V + LR_apply(Fs2r,U_cell{2});
    Utmp3 = HSS1D_apply_fsym(FI.rs_Inv,Utmp3); 
    U_cell{3} = LR_apply(T,Utmp3); 
    
    U3d = LR_apply(T,FId_rs_Inv*(T.V+LR_apply(Fs2r,V1d))); 
    
    % U*V.' = U1*V1.'+U2*V2.'+U3*V3.'
    [U,V] = LR_recomp(U_cell,V_cell,acc/10);
    
    % HSS1D compression of lowrank 
    M = HSS1D_compressLR_fsym(U,V,FI.SrsInv,acc);
    M = HSS1D_symmetrize_fsym(M); 
    
    Md = U1d*V1d.'+V1d*U1d.'+U3d*V3d.'; 
    
    fprintf('\n Low Rank Perturbation error \n')
    E_LR = norm(Md*Vsk - HSS1D_apply_fsym(M,Vsk))./norm(Md*Vsk); 
    display(E_LR)
    
    % EI = SrsInv + M
    EI = HSS1D_sum_fsym(FI.SrsInv,M,acc); 
    EI = HSS1D_recompress_fsym(EI,acc);
    EI = HSS1D_fsym_to_nsym(EI); 
    
    EId = FId_SrsInv + Md;
    %Ed = inv(EId); 
    
    % E = inv(EI)
    [E_HSSinv,E_FTinv] = HSS1D_invert_nsym(EI,acc);
    E = HSS1D_transforminv_nsym(EI,E_HSSinv,E_FTinv);       
    E = HSS1D_nsym_to_fsym(E,acc); 
    
    fprintf('\n E matrix test error \n')
    E_inv = norm(EId\Vsk - HSS1D_apply_fsym(E,Vsk))./norm(EId\Vsk); 
    display(E_inv)
    
else

    % (I) EInv = phi_sk + T*phi_rs2sk + phi_sk2rs*T + T*phi_rs*R_dn,
    % where phis are defined by the block inversion formulae. 
    if (opts.LR==0)
        phi_rs = FI.rs_Inv*(eye(length(FI.rs_Inv)) + FI.sk2rs*FI.SrsInv*FI.rs2sk*FI.rs_Inv);
        phi_rs2sk = -FI.SrsInv*FI.rs2sk*FI.rs_Inv; 
        phi_sk2rs = -FI.rs_Inv*FI.sk2rs*FI.SrsInv; 
        
        EI = FI.SrsInv + T*phi_sk2rs + phi_rs2sk*T.' + T*phi_rs*T.'; 
        
        % (II) Inversion of EInv, direct
        %E = inv(EI); 
        E = LOCAL_Invert_sym(EI); 
    else
        
        J_F = FI.J; 
        Ind = 1:length(J); Jinv(J) = Ind; 
        
        J_J2F = Jinv(J_F);
        Js_J2F = J_J2F(1:k); 
        Jr_J2F = J_J2F(k+1:end) - k; 
        
        % Permute T to make indices the same as in F
        T.U = T.U(Js_J2F,:); 
        T.V = T.V(Jr_J2F,:); 
        
        U_cell = cell(3,1); V_cell = cell(3,1); 
        V_cell{1} = FI.SrsInv*LR_apply(FI.rs2sk,FI.rs_Inv*(T.V)); 
        U_cell{1} = -T.U;
        
        V_cell{2} = U_cell{1};          
        U_cell{2} = V_cell{1}; 
        
        V_cell{3} = T.U;          
        Fs2r.U = FI.rs2sk.V; Fs2r.V = FI.rs2sk.U; 
        U_cell{3} = T.U*T.V.'*(FI.rs_Inv*(T.V + LR_apply(Fs2r,U_cell{2}))); 
        
        [U,V] = LR_recomp(U_cell,V_cell,acc); 
        
        EI = FI.SrsInv + U*V.';
        E = LOCAL_Invert_sym(EI); 
        
        % Inversion using SMW
        %S = FI.SrsInv\eye(k); 
        %m = size(U,2); 
        %MI = (eye(m) + V'*(FI.SrsInv\U)); 
        %E = S*(eye(k) - U*(MI\V')*S); 
    end
end
end

function Ai = LOCAL_Invert_sym(A)
 
Ai = inv(A);
Ai = (Ai+Ai.')/2; 

end