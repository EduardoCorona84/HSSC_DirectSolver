function [E,Js_J2F] = HSS2D_build_E_nsym(FI,J,T_up,T_dn,k,acc,opts)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%     FUNCTION CALL:
%         [E,Js_J2F] = HSS2D_build_E_nsym(FI,J,T_up,T_dn,k,acc,opts)
%
%      DESCRIPTION:
%         This function is an implementation of Algorithm 5 in the paper - which corresponds to
%         line 10 of Algorithm 2. For a given box, it computes E as a function of the
%         interpolation operators and the blocks of F^-1. It makes no assumptions on the
%         symmetry of the associated kernel (nsym).
%
%     INPUT:
%         FI      <nxn or struct>     (output of HSS2D_build_Finv_fsym.m) The blocks of F^-1
%                                         needed to compute E.
%         J       <kx1 int>           [I_sk I_rs] Permutation vector that places skeleton points
%                                         first, residual points second.
%         T       <kxk or lr>         Interpolation matrix (stored either densely or as a
%                                         rank (U,V) pair).
%         k       <int>               Size of skeleton.
%         acc     <float>             Desired accuracy.
%
%         opts    <struct>            .mat determines how E is stored. If 'dense' or 'block', E
%                                         is dense and .mat.LR specifies if the interpolation
%                                         operations are stored as low-rank. If .mat = 'HSS', E
%                                         is HSS1D and interpolation operators are assumed to be
%                                         low-rank.
%
%     OUTPUT:
%         E   <kxk or HSS1D>   Schur Complement matrix used to represent A^-1.
%
%                         EInv = R*FInv*L = [I T] * FI * [I ; T.']
%                         E = inv(EInv)
%




if strcmp(opts.mat,'dense')
    
    if (opts.LR == 0)
        R(:,J) = [ eye(k) T_up.U]; 
        L(J,:) = [ eye(k) T_dn.U].'; 
        EI = R*(FI*L);
        
        % Direct inversion
        E = inv(EI); 
    else
        Phi = FI(J,J);
        U_cell = cell(3,1); V_cell = cell(3,1); 
        V_cell{1} = Phi(k+1:end,1:k).'*T_up.V;  U_cell{1} = T_up.U;
        V_cell{2} = T_dn.U;                    U_cell{2} = Phi(1:k,k+1:end)*(T_dn.V); 
        V_cell{3} = T_dn.U;                    
        
        U_cell{3} = T_up.U*(T_up.V.')*Phi(k+1:end,k+1:end)*T_dn.V; 
        %[U,V] = LR_recomp(U_cell,V_cell,acc); 
        
        % Inversion
        EI = Phi(1:k,1:k) + U_cell{1}*V_cell{1}.' + U_cell{2}*V_cell{2}.' + U_cell{3}*V_cell{3}.'; 
        E = inv(EI); 
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
        T_up.U = T_up.U(Js_J2F,:); T_dn.U = T_dn.U(Js_J2F,:); 
        T_up.V = T_up.V(Jr_J2F,:); T_dn.V = T_dn.V(Jr_J2F,:); 
        
        U_cell = cell(3,1); V_cell = cell(3,1); 
        
        % U1,V1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % V1.' = T_up.V.' * FrsInv * F_s2r * SrsInv, U1 = -T_up.U; 
        
        U_cell{1} = -T_up.U;          
        
        FT = HSS1D_transpose(FI.rs_Inv); 
        ST = HSS1D_transpose(FI.SrsInv); 
        Vtmp1 = HSS1D_apply_nsym(FT,T_up.V); 
        Ft_s2r.U = FI.sk2rs.V; Ft_s2r.V = FI.sk2rs.U; 
        Vtmp1 = LR_apply(Ft_s2r,Vtmp1); 
        V_cell{1} = HSS1D_apply_nsym(ST,Vtmp1); 
        
        % U2,V2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % V2 = -T_dn.U, U2 = SrsInv * F_r2s * FrsInv * T_dn.V ; 
        Utmp2 = HSS1D_apply_nsym(FI.rs_Inv,T_dn.V); 
        Utmp2 = LR_apply(FI.rs2sk,Utmp2); 
        
        U_cell{2} = HSS1D_apply_nsym(FI.SrsInv,Utmp2); 
        V_cell{2} = -T_dn.U;
        
        % U3,V3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % U3 = T_up * Frs_Inv * (T_dn.V + F_s2r * V2), V3 = T_dn.U; 
        V_cell{3} = T_dn.U;  
        
        Utmp3 = T_dn.V + LR_apply(FI.sk2rs,U_cell{2});
        Utmp3 = HSS1D_apply_nsym(FI.rs_Inv,Utmp3); 
        U_cell{3} = LR_apply(T_up,Utmp3); 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % U*V.' = U1*V1.'+ U2*V2.' + U3*V3.'
        [U,V] = LR_recomp(U_cell,V_cell,acc/10);
        M = HSS1D_compressLR_nsym(U,V,FI.SrsInv,acc);
        % EI = SrsInv + U*V.'
        EI = HSS1D_sum_nsym(FI.SrsInv,M,acc); 
        EI = HSS1D_recompress_nsym(EI,acc);
        
        [E_HSSinv,E_FTinv] = HSS1D_invert_nsym(EI,acc);
        E = HSS1D_transforminv_nsym(EI,E_HSSinv,E_FTinv);
    
else

    % (I) EInv = phi_sk + T*phi_rs2sk + phi_sk2rs*T + T*phi_rs*R_dn,
    % where phis are defined by the block inversion formulae. 
    if (opts.LR==0)
        phi_rs = FI.rs_Inv*(eye(length(FI.rs_Inv)) + FI.sk2rs*FI.SrsInv*FI.rs2sk*FI.rs_Inv);
        phi_rs2sk = -FI.SrsInv*FI.rs2sk*FI.rs_Inv; 
        phi_sk2rs = -FI.rs_Inv*FI.sk2rs*FI.SrsInv; 
        
        EI = FI.SrsInv + T_up*phi_sk2rs + phi_rs2sk*T_dn.' + T_up*phi_rs*T_dn.'; 
        
        % (II) Inversion of EInv, direct
        E = inv(EI); 
    else
        
        J_F = FI.J; 
        Ind = 1:length(J); Jinv(J) = Ind; 
        
        J_J2F = Jinv(J_F);
        Js_J2F = J_J2F(1:k); 
        Jr_J2F = J_J2F(k+1:end) - k; 
        
        % Permute T to make indices the same as in F
        T_up.U = T_up.U(Js_J2F,:); T_dn.U = T_dn.U(Js_J2F,:); 
        T_up.V = T_up.V(Jr_J2F,:); T_dn.V = T_dn.V(Jr_J2F,:); 
        
        U_cell = cell(3,1); V_cell = cell(3,1); 
        % U1,V1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % V1.' = T_up.V.' * FrsInv * F_s2r * SrsInv, U1 = -T_up.U; 
        V_cell{1} = (FI.SrsInv.'*LR_apply(FI.rs2sk,FI.rs_Inv.'*T_up.V)).';          
        U_cell{1} = -T_up.U;
        % U2,V2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % V2 = -T_dn.U, U2 = SrsInv * F_r2s * FrsInv * T_dn.V ; 
        U_cell{2} = FI.SrsInv*LR_apply(FI.rs2sk,FI.rs_Inv*(T_dn.V)); 
        V_cell{2} = -T_dn.U;
        % U3,V3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % U3 = T_up * Frs_Inv * (T_dn.V + F_s2r * V2), V3 = T_dn.U; 
        V_cell{3} = T.U;          
        U_cell{3} = T.U*T.V.'*(FI.rs_Inv*(T.V + LR_apply(FI.sk2rs,U_cell{2}))); 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % U*V.' = U1*V1.'+ U2*V2.' + U3*V3.'
        [U,V] = LR_recomp(U_cell,V_cell,acc); 
        
        EI = FI.SrsInv + U*V.';
        E = inv(EI); 
    end
end
end