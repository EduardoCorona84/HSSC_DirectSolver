function Y = HSS2D_applyFinv_fsym(FI,V,k,opts)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%     FUNCTION CALL:
%         Y = HSS2D_applyFinv_fsym(FI,V,k,opts)
%
%     DESCRIPTION:
%         This function applies F^-1 to a vector, where F is the 2x2 compressed block matrix as
%         defined in the paper. The Schur block inversion formulae for the 4 blocks of F^-1 are
%         applied to V(FI.J), and then we invert the permutation to get Y. This function assumes
%         the underlying kernel matrix is full-symmetric (fsym).
%
%
%     INPUT:
%         FI      <struct>    (output of HSS2D_build_Finv_fsym.m) contains blocks needed to
%                                 apply Finv. FI.J holds the indices of skeleton followed by
%                                 residual points.
%         V       <nx1 float> Column vector defined on I_src.
%         k       <int>       Size of skeleton I_sk.
%         opts    <string>    {'dense', 'block', 'HSS'} Specifies how E,F^-1 are stored. It is
%                                 determined by the level in the tree, params.lev_HSS, and
%                                 params.Fopts in the parent function HSS2D_applyinv_fsym.
%                 - 'dense'   E,F^-1 are each computed as a single dense matrix.
%                 - 'block'   The blocks of FI.* are computed densely. This corresponds to the
%                                 left side of Algorithm 4.
%                 - 'HSS'     The blocks of F^-1 are HSS1D or low-rank. E matrices are HSS1D.
%                                 This corresponds to the right side of Algorithm 4.
%
%
%     OUTPUT:
%         Y       <nx1 float>     = FI*V  defined on I_src
%




%Y = zeros(size(J_up,1),1); 

if strcmp(opts,'dense')
    % Dense apply
    Y = FI*V; 
    
elseif strcmp(opts,'HSS')
    % HSS / LR block apply
    J = FI.J; 
    
    %permute V so that it is [V(I_sk) V(I_rs)]
    V = V(J,:); 

    V_sk = V(1:k,:);
    V_rs = V(k+1:end,:); 

    %We now compute the 4 products needed using the Schur complement formulae

    % (I) Result for I_sk: 
    r = length(J)-k; 

    if r>0
        % HSS apply
        psi = HSS1D_apply_fsym(FI.rs_Inv,V_rs); 
        % LR apply
        eta = V_sk - LR_apply(FI.rs2sk,psi); 
    else
        eta = V_sk;
    end

    % Y defined at skeleton points I_sk
    %HSS apply
    Y_sk = HSS1D_apply_fsym(FI.SrsInv,eta); 

    % (II) Result for I_rs
    if r>0
        % LR apply
        F_s2r.U = FI.rs2sk.V; F_s2r.V = FI.rs2sk.U; 
        lam =  LR_apply(F_s2r,Y_sk);   

        %HSS apply
        Y_rs = HSS1D_apply_fsym(FI.rs_Inv,(V_rs - lam)); 

        % (III) form and permute Y.
        Y = [Y_sk ; Y_rs];
    else
        Y = Y_sk;  
    end
    
    Y(J,:) = Y;  
    
    
else
    % Dense block apply
    J = FI.J;
    
    %permute V so that it is [V(I_sk) V(I_rs)]
    V = V(J,:); 

    V_sk = V(1:k,:);
    V_rs = V(k+1:end,:); 

    %We now compute the 4 products needed using the Schur complement formulae

    % (I) Result for I_sk: 

    % HSS apply
    psi = FI.rs_Inv*V_rs; 
    % LR apply
    eta = LR_apply(FI.rs2sk,psi); 

    % Y defined at skeleton points I_sk
    %HSS apply
    Y_sk = FI.SrsInv*(V_sk - eta); 

    % (II) Result for I_rs
    % LR apply
    F_s2r.U = FI.rs2sk.V; F_s2r.V = FI.rs2sk.U; 
    lam =  LR_apply(F_s2r,Y_sk);   

    %HSS apply
    Y_rs = FI.rs_Inv*(V_rs - lam); 

    % (III) form and permute Y.
    Y = [Y_sk ; Y_rs];
    Y(J,:) = Y; 
    
end
%--------------------------------------------------------------------------
