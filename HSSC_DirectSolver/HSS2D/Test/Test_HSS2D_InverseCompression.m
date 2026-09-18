function [Err_DvsC,Err_D,Err_C,T_tree,TD_inv,TC_inv,TD_apply,TC_apply,MB_AD,MB_AC] = Test_HSS2D_InverseCompression(params,INFO_TREE)
% Tests Inverse Compression and Apply for a single instance, determined by the 
% parameters on the params struct. If the corresponding binary tree 
% INFO_TREE is included as a second input, tree build is skipped and 
% INFO_TREE is used for inverse compression. 

% For the apply, we Err_Dsure relative error of applying the compressed 
% HSS2D inverses AD_HSS,AC_HSS to randomly generated vectors.

%{
 Outputs:
 INFO_TREE - HSS2D binary tree
 Err_D     - mean error for dense block algorithm      (HSS-D)
 Err_C     - mean error for compressed block algorithm (HSS-C)
 Err_DvsC  - mean error between HSS-D and HSS-C
 T_tree    - time to build binary tree
 TD_comp   - time to compress HSS-D inverse AD_HSS
 TC_comp   - time to compress HSS-C inverse AC_HSS
 TD_apply  - mean apply time (HSS-D)
 TC_apply  - mean apply time (HSS-C)
%}

Test_setup; 

rng('default');     % For older versions of MATLAB, use line below instead
%rand( 'seed',0);   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ntot = size(params.X_source,1); % Ng^2*np^2;
TI = params.transinv; 
sym = params.sym; 
% Maximum N for which dense-block algorithm is run 
% (if there are memory constraints)
Ndense = 1000000; 

% If there's only one argument (params), build HSS tree. 
if nargin == 1
    fprintf('\n Parameters \n'); 
    display(params)  

    % 2D binary tree build & skeleton info
    fprintf('\n (0) HSS2D binary tree build & skeleton info \n'); 
    if params.sym == 0
        tic; 
        [INFO_TREE]=HSS2D_build_tree_nsym(params);
        T_tree = toc; 
    else
        tic; 
        [INFO_TREE]=HSS2D_build_tree_fsym(params); 
        T_tree = toc; 
    end

    fprintf('\n Tree Build Time: %e \n',T_tree);
else
    T_tree = 0; 
end

% Determine lev_HSS (level to switch from dense to compressed blocks)
if length(INFO_TREE.BOX)>7
    n_cut = params.n_cut; 
    lev_HSS = levHSS(INFO_TREE,n_cut,params.transinv);
    lev_HSS = min(lev_HSS,INFO_TREE.depth-2); 
else
    lev_HSS = -1; 
end

fprintf('\n First level of compressed blocks lev_HSS: \n'); 
display(lev_HSS)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverse Compression
fprintf('\n (I) HSS2D Inverse Compression \n'); 

X = params.X_source;

if params.dim == 3
    fprintf('\n 3D shape of type: \n')
    display(params.type)  
    params.T_source = params.X_source;    
   params.X_source = shape_3D(params.T_source,params,params.type);   
   display(params.type)
end

if params.sym == 0
    fprintf('\n (I.1) Dense - O(N^1.5) \n'); 
    params.Fopts = 'dense';
    
    if Ntot<Ndense
    tic; 
    [A_HSSD] = HSS2D_invert_nsym(INFO_TREE,params);
    TD_inv = toc; 
    fprintf('Inversion Time: %e \n',TD_inv)
    
    TD_inv = mean(TD_inv); 
    else
        TD_inv = 0; 
    end
    
    fprintf('\n (I.2) HSS - O(N) \n'); 
    params.Fopts = 'HSS'; params.lev_HSS = lev_HSS; 
    
    tic; 
    [A_HSSC] = HSS2D_invert_nsym(INFO_TREE,params);
    TC_inv = toc; 
    fprintf('Inversion Time: %e \n',TC_inv)    
    
    TC_inv = mean(TC_inv); 
    
else
    if Ntot<Ndense
        fprintf('\n (I.1) Dense - O(N^1.5) \n');
        params.Fopts = 'dense';
        
        tic; 
        [A_HSSD] = HSS2D_invert_fsym(INFO_TREE,params);
        TD_inv = toc; 
        fprintf('\n HSS-D Inversion Time: %e \n',TD_inv)
    else
       TD_inv = 0;   
    end
    
    TD_inv = mean(TD_inv); 
    
    fprintf('\n (I.3) HSS - O(N) \n'); 
    params.Fopts = 'HSS'; params.lev_HSS = lev_HSS;                   
    display(params.Fopts);    
    
    tic; 
    [A_HSSC] = HSS2D_invert_fsym(INFO_TREE,params);
    TC_inv = toc; 
    fprintf('\n HSS-C Inversion Time: %e \n',TC_inv)
    
    TC_inv = mean(TC_inv); 
end

fprintf('\n (III) Inverse Apply \n'); 
Ntot = size(X,1); 
m = 1; 
apply_params.transinv = params.transinv; 
apply_params.INTERform = params.INTERform; 

if params.kh == 0
    fprintf('\n Non-oscillatory problem, random rhs \n'); 
    Vec = randn(Ntot,m); 
    Vec = Vec.*(ones(Ntot,1)*(1./sqrt(sum(Vec.*Vec))));
else
    fprintf('\n Lippman Schwinger, scattered incoming field \n'); 
    params.Fopts = 'dense'; apply_params.Fopts = 'dense'; 
    Vec = exp(1i*params.kh*X(:,1)); 
    Bump = Bump_function(X(:,1),X(:,2),4); 
    
    if sym == 0
        M_HSS = HSS2D_compress_nsym(INFO_TREE,params); 
        Vec = -Bump.*Vec; 
    else
        M_HSS = HSS2D_compress_fsym(INFO_TREE,params); 
        Vec = -sqrt(Bump).*Vec; 
    end
    %{
    w = (1 + params.kh^2*(0.25*1i)*params.dr_weights);
    
    if sym == 0
        M_HSS = HSS2D_compress_nsym(INFO_TREE,params); 
        Vec = -(HSS2D_apply_nsym(INFO_TREE,M_HSS,Bump.*Vec,apply_params) - w*Bump.*Vec); 
    else
        M_HSS = HSS2D_compress_fsym(INFO_TREE,params); 
        Vec = -(HSS2D_apply_fsym(INFO_TREE,M_HSS,Bump.*Vec,apply_params) - w*Bump.*Vec); 
    end
    %}
end

if params.sym == 0
    % Direct Solution
    fprintf('\n Direct Solution \n')
    apply_params.Fopts = 'dense';
    if Ntot<Ndense
        tic; 
        Y_dense = HSS2D_applyinv_nsym(INFO_TREE,A_HSSD,Vec,apply_params);
        TD_apply = toc/m;
    else
        TD_apply = 0; 
    end
    display(TD_apply)
    
    apply_params.Fopts = 'HSS'; apply_params.lev_HSS = lev_HSS; 
    tic; 
    Y_HSS   = HSS2D_applyinv_nsym(INFO_TREE,A_HSSC,Vec,apply_params);
    TC_apply = toc/m; 
    display(TC_apply)
    
    if params.kh == 0
        params.Fopts = 'dense';
        if sym == 0
            M_HSS = HSS2D_compress_nsym(INFO_TREE,params); 
        else
            M_HSS = HSS2D_compress_fsym(INFO_TREE,params);  
        end
    end
    
    % Preconditioned Iterative Solution
    fprintf('\n Iterative Solution \n')
    paramsfs = params; paramsfs.sym = 1; paramsfs.transinv = 1; paramsfs.layers = 2; 
    apply_paramsfs = apply_params; apply_paramsfs.transinv = 1; 
    
    display(paramsfs)
    
    tic; 
    TREEfs =HSS2D_build_tree_fsym(paramsfs); 
    paramsfs.lev_HSS = levHSS(TREEfs,paramsfs.n_cut,1); paramsfs.Fopts = 'HSS';
    [I_HSS] = HSS2D_invert_fsym(TREEfs,paramsfs);
    Tinv_iter = toc; 
    
    display(Tinv_iter)
    
    apply_params.Fopts = 'dense'; 
    Mapp = @(x) HSS2D_apply_nsym(INFO_TREE,M_HSS,x,apply_params); 
    apply_paramsfs.Fopts = 'HSS'; apply_paramsfs.lev_HSS = paramsfs.lev_HSS;
    display(apply_paramsfs)
    Iapp = @(x) HSS2D_applyinv_TI(TREEfs,I_HSS,x,apply_paramsfs); 
    
    tic; 
    [Y_iter,Flag,Relres,IT] = bicgstab(Mapp,Vec,params.acc/10,200,Iapp); 
    Titer_apply = toc; 
    
    fprintf('Number of iterations: %d \n',IT); 
    display(Titer_apply)
    display(Relres)
    
else
    apply_params.Fopts = 'dense';
    if Ntot<Ndense
        tic; 
        if TI == 1
        Y_dense = HSS2D_applyinv_TI(INFO_TREE,A_HSSD,Vec,apply_params);
        else
        Y_dense = HSS2D_applyinv_fsym(INFO_TREE,A_HSSD,Vec,apply_params);    
        end
        TD_apply = toc/m; 
    else
        TD_apply = 0; 
    end
        
    apply_params.Fopts = 'HSS'; apply_params.lev_HSS = lev_HSS;    

    tic; 
    if TI == 1
       Y_HSS   = HSS2D_applyinv_TI(INFO_TREE,A_HSSC,Vec,apply_params);
    else
       Y_HSS   = HSS2D_applyinv_fsym(INFO_TREE,A_HSSC,Vec,apply_params); 
    end
    TC_apply = toc/m;
end

     
    params.Fopts = 'dense'; 
    if params.sym == 1
        [M_HSS] = HSS2D_compress_fsym(INFO_TREE,params);
    else
        %[M_HSS] = HSS2D_compress_nsym(INFO_TREE,params);
    end
    
    apply_params.Fopts = 'dense'; 
    if TI==1
        if Ntot<Ndense
        Vd      = HSS2D_apply_TI(INFO_TREE,M_HSS,Y_dense,apply_params);
        end
        VH      = HSS2D_apply_TI(INFO_TREE,M_HSS,Y_HSS,apply_params);
    elseif sym == 1
        if Ntot<Ndense
        Vd      = HSS2D_apply_fsym(INFO_TREE,M_HSS,Y_dense,apply_params);
        end
        VH      = HSS2D_apply_fsym(INFO_TREE,M_HSS,Y_HSS,apply_params);
    else
        if Ntot<Ndense
        Vd      = HSS2D_apply_nsym(INFO_TREE,M_HSS,Y_dense,apply_params);
        end
        VH      = HSS2D_apply_nsym(INFO_TREE,M_HSS,Y_HSS,apply_params);
    end
    
    % Iterative refinement (test)
    R = Vec - VH; 
    
    apply_params.Fopts = 'HSS'; 
    if TI == 1
        Y2_HSS = Y_HSS + HSS2D_applyinv_TI(INFO_TREE,A_HSSC,R,apply_params);
    elseif params.sym == 1 
        Y2_HSS = Y_HSS + HSS2D_applyinv_fsym(INFO_TREE,A_HSSC,R,apply_params);    
    else
        Y2_HSS = Y_HSS + HSS2D_applyinv_nsym(INFO_TREE,A_HSSC,R,apply_params);    
    end
    
    apply_params.Fopts = 'dense'; 
    if params.sym == 1
        V2H      = HSS2D_apply_fsym(INFO_TREE,M_HSS,Y2_HSS,apply_params);
    else
        V2H      = HSS2D_apply_nsym(INFO_TREE,M_HSS,Y2_HSS,apply_params);
    end

if Ntot<10000   
    if params.dim == 3
        fprintf('\n 3D: K(X,X) (dense) vs HSS2D \n');     
        X = shape_3D(X,params,params.type);
        display(params.type)
    end
    A = Kernel_Eval(X,X,params);
    Y_true = A\Vec;
    Errd_Apply = sqrt(max(sum((Y_true-Y_dense).^2)));
    ErrH_Apply  = sqrt(max(sum((Y_true-Y_HSS).^2)));
    Err_HIR = sqrt(max(sum((Y_true-Y2_HSS).^2)));
    Err_HSS    = sqrt(max(sum((Y_dense-Y_HSS).^2))); 
    %Errb_Apply  = sqrt(max(sum((Y_true - Y_block).^2))); 
elseif Ntot<10000000
    
    if Ntot<Ndense
        Errd_Apply = sqrt(max(sum((Vec-Vd).^2)));
        Err_HSS    = norm(Y_dense-Y_HSS)/norm(Y_dense);
    else
        Errd_Apply = 0; 
        Err_HSS = 0;
    end
    
    ErrH_Apply = sqrt(max(sum((Vec-VH).^2)));
    Err_HIR = sqrt(max(sum((Vec-V2H).^2)));
    
    %Errb_Apply = sqrt(max(sum((Vec-VB).^2))); 
else
    Errd_Apply = 0;
    ErrH_Apply = 0;
    Err_HSS    = 0;
    Err_HIR = 0; 
    
end

Err_DvsC = Err_HSS; 
Err_D = Errd_Apply;  
Err_C = ErrH_Apply; 

if Ntot<Ndense
info_AD = whos('A_HSSD'); MB_AD = info_AD.bytes/(2^20); RpDOF_AD = info_AD.bytes/(8*Ntot); 
end
info_T  = whos('INFO_TREE'); MB_T = info_T.bytes/(2^20); RpDOF_T = info_T.bytes/(8*Ntot); 
info_AC = whos('A_HSSC'); MB_AC = info_AC.bytes/(2^20); RpDOF_AC = info_AC.bytes/(8*Ntot); 

fprintf('\n mean Error of HSS vs Dense Inverse Application: %e \n',Err_DvsC);
fprintf('\n mean Error of Dense Inverse Application: %e \n',Err_D);
fprintf('\n mean Error of HSS Inverse Application: %e \n',Err_C);
fprintf('\n mean Error of HSS Inverse Application (after 1 IR): %e \n',Err_HIR);
fprintf('\n mean O(N^1.5) Apply Time: %e \n',TD_apply);
fprintf('\n mean O(N) Apply Time: %e \n',TC_apply);


fprintf(1,'\n Memory required for O(N) Tree build = %7.2f MB     (about %0d reals/DOF) \n',...
          MB_T,RpDOF_T);   
if Ntot<Ndense
fprintf(1,'\n Memory required for O(N^1.5) Inverse = %7.2f MB     (about %0d reals/DOF) \n',...
          MB_AD,RpDOF_AD);
end
fprintf(1,'\n Memory required for O(N) Inverse = %7.2f MB     (about %0d reals/DOF) \n',...
          MB_AC,RpDOF_AC);   
      
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%