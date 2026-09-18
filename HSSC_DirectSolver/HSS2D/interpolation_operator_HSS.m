function T = interpolation_operator_HSS(cB,lev,par,X_proxy,X1,X2,type,q)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%     FUNCTION CALL:
%          T = interpolation_operator_HSS(cB,lev,par,X_proxy,X1,X2,type,q)
%
%     DESCRIPTION:
%         This function calculates the interpolation operators T that comprise L and R in the
%         telescoping factorization. To speedup computation we use equivalent densities (proxy
%         set), we compress K(X_proxy,X1) as HSS1D, and we compress K(X_proxy,X2) as low-rank.
%         We solve for T by HSS1D least squares:
%
%                     T(X1,X2) = K(X_proxy,X1) \ K(X_proxy,X2)
%
%         In other words, for
%
%                     phi_1 = T(X1,X2) phi_2
%                     K(X_proxy,X1) phi_1 = K(X_proxy,X2) phi_2
%
%         This function corresponds to Algorithm 3 in the paper. It differs from
%         interpolation_operator.m by using HSS1D compression and least-squares. It is called
%         within HSS2D_bintree_* during the BUILD TREE stage for boxes with k > par.n_cut.
%
%     INPUT:
%         cB      <2x1 float>         Center of the box.
%         lev     <int>               Level in the tree.
%         par     <struct>            (see HSS_tree_parameters.m) Solve parameters.
%         X_proxy <mx2 float>         Proxy points used in place of X_ext.
%         X1      <kx2 float>         Skeleton points.
%         X2,     <(n-k)x2 float>     Residual points.
%         type    <string>            {'box_2_proxy', 'proxy_2_box'} Direction of interpolation.
%         q       <(depth+1)x1 int>   Array of interpolation operator ranks (used in randomized
%                                     ID to produce a first guess).
%
%         OUTPUT:
%             T[X1,X2] <kxq,qx(n-k) float>    Interpolation operator stored as low-rank UV'
%                                             pair, where q is the rank.
%




nbox_max = 108; 
acc = par.acc; 
tau = par.tau; 
mu = par.mu; 
maxit = 3;     

% Check to see if the inputs are non-empty
if (size(X1,1) > 0 && size(X2,1) > 0)
    
    % If skeleton is not previously sorted, it needs to be sorted around a
    % curve. 
    %[X_bds,Jbd,~,~,~] = sort_skeleton(X1,cB,lev,h,layers);
    
    %T is then K(X_proxy,X1)^(-1) * K(X_proxy,X2)
    %We compress K21 and use this to speed up the solve (since it is of low rank)
    if strcmp(type,'box_2_proxy')
        % K[X_ext,X1], HSS compression of proxy-skeleton interactions
        KH11 = HSS1D_compress_rectangular_nsym('green',X_proxy,X1,nbox_max,acc/10,par);
        
        % Alternative if matrix is square: 
        %{
        KH11 = HSS1D_compress_green_nsym_Kproxy(X_proxy,X1,nbox_max,1e-11,par);
        [HSSinv,FTinv] = HSS1D_invert_nsym(KH11,1e-11);
    
        % Function handles for matrix and inverse apply
        mvK11 = @(x) HSS1D_apply_nsym(KH11,x,1); 
        mvK11I = @(x) HSS1D_applyinv_nsym(KH11,HSSinv,FTinv,x,1); 
        %}
        
        % K_[X_ext,X2]
        K21 = Kernel_Eval(X_proxy,X2,par); 
        
        % Low Rank decomposition using Randomized ID and dense
        % matvecs: 
        mvK21 = @(x) K21*x; mvK21t = @(x) K21.'*x;
        m = size(K21,1); n = size(K21,2); 
        
        % Randomized ID parameters. q is an array storing ranks of
        % interpolation operators in finer levels. See interpolation
        % operator for details about rank estimation. 
        if q(lev+1)>0
            ID_params.q = q(lev+1)+1; 
            ID_params.C = 10;
        else
            ID_params.q = min(ceil(2*q(lev+3) - q(lev+5) + 2),min(m,n)); 
            ID_params.C = ceil(q(lev+3) - q(lev+5)/2) + 1;
        end
        
        ID_params.par = acc; ID_params.form = 'dense'; 
        
        % Randomized ID
        K21_LR = LRID_rand(m,n,mvK21,mvK21t,'norm',ID_params);
               
        % T = K11\K21 (HSS least squares)
        T.U = HSS1D_OLS(KH11,K21_LR.U,acc,tau,mu,maxit,'qr');
        T.V = K21_LR.V;
        
        % Alternative to least squares if matrix is square: a few rounds of
        % bicgstab or iterative refinement to improve accuracy of
        % ill-conditioned solve. 
        %{
        [T.U,~] = bicgstabvec(mvK11,K21_LR.U,acc,10,mvK11I);  
        %T.U = Iterative_Refinement(mvK11I,mvK11,K21_LR.U,acc,10); 
        T.U(Jbd,:) = T.U; 
        %}
        
    else
        % K[X_ext,X1], HSS compression of skeleton-proxy interactions
        KH11 = HSS1D_compress_rectangular_nsym('green',X1,X_proxy,nbox_max,acc/10,par);
        KH11 = HSS1D_transpose_rectangular(KH11); 
        
        % Alternative if matrix is square: 
        %{
        KH11 = HSS1D_compress_green_nsym_Kproxy(X1,X_proxy,nbox_max,1e-11,par);
        KH11 = HSS1D_transpose(KH11); 
        [HSSinv,FTinv] = HSS1D_invert_nsym(KH11,1e-11);
    
        % Function handles for matrix and inverse apply
        mvK11 = @(x) HSS1D_apply_nsym(KH11,x,1); 
        mvK11I = @(x) HSS1D_applyinv_nsym(KH11,HSSinv,FTinv,x,1); 
        %}
        
        % K_[X2,X_ext]
        K21 = Kernel_Eval(X2,X_proxy,par); 
        
        % Low Rank decomposition using Randomized ID and dense
        % matvecs: 
        mvK12 = @(x) K21.'*x; mvK12t = @(x) K21*x;
        m = size(K21',1); n = size(K21',2); 
        
        % Randomized ID parameters. q is an array storing ranks of
        % interpolation operators in finer levels. 
        if q(lev+1)>0
            ID_params.q = q(lev+1)+1; 
            ID_params.C = 10;
        else
            ID_params.q = min(ceil(2*q(lev+3) - q(lev+5) + 2),min(m,n)); 
            ID_params.C = ceil(q(lev+3) - q(lev+5)/2) + 1;
        end
        
        ID_params.par = acc; ID_params.form = 'dense'; 
        
        % Randomized ID
        K12_LR = LRID_rand(m,n,mvK12,mvK12t,'norm',ID_params);
        
        % T = K11\K12 (HSS least squares)
        T.U = HSS1D_OLS(KH11,K12_LR.U,acc,tau,mu,maxit,'qr');   
        T.V = K12_LR.V;
        
         % Alternative to least squares if matrix is square: a few rounds of
        % bicgstab or iterative refinement to improve accuracy of
        % ill-conditioned, HSS square matrix solve. 
        %{
        [T.U,~] = bicgstabvec(mvK11,K12_LR.U,acc,10,mvK11I);  
        %T.U = Iterative_Refinement(mvK11I,mvK11,K12_LR.U,acc,10); 
        T.U(Jbd,:) = T.U; 
        %}
        
    end
    

else
    
    % Returns empty array
    T = []; 

end

end