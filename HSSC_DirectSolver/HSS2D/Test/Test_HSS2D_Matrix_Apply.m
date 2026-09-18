function [Err_D,Err_C,Err_DvsC,T_tree,TD_comp,TC_comp,TDA,TCA] = Test_HSS2D_Matrix_Apply(params,INFO_TREE)
% Tests Matrix Compression and Apply for a single instance, determined by the 
% parameters on the params struct. If the corresponding binary tree 
% INFO_TREE is included as a second input, tree build is skipped and 
% INFO_TREE is used for matrix compression. 

% For the apply, we measure relative error of applying the compressed 
% HSS2D matrices MD_HSS,MC_HSS to randomly generated vectors.

%{
 Outputs:
 INFO_TREE - HSS2D binary tree
 Err_D       - mean error for dense block algorithm      (HSS-D)
 Err_DvsC       - mean error for compressed block algorithm (HSS-C)
 T_tree    - time to build binary tree
 T_comp    - time to compress matrix A_HSS
 TDA      - mean apply time (HSS-D)
 TCA      - mean apply time (HSS-C)
%}

Test_setup; 

rng('default');     % For older versions of MATLAB, use line below instead
%rand( 'seed',0); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%params = HSS_tree_parameters(pot,kh,sym,Ng,np,acc,lay,n_cut,TI);
%Ntot = Ng^2*np^2;
Ntot = size(params.X_source,1); 
TI = params.transinv; 

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

fprintf('\n Tree Build Time O(N): %e \n',T_tree);

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
% Matrix Compression
fprintf('\n (I) HSS2D Matrix Compression \n'); 

X = params.X_source;

if params.dim == 3
    fprintf('\n 3D shape of type: \n')
    display(params.type)  
   params.T_source = params.X_source;    
   params.X_source = shape_3D(params.T_source,params,params.type);    
   display(params.type)
end

% Matrix Compression using dense blocks (MD_HSS) and compressed blocks
% (MC_HSS).
if params.sym == 0    
    if Ntot<Ndense
        fprintf('\n (I.1) Dense - O(Nlog(N)) \n'); 
        params.Fopts = 'dense';
    
        tic; 
        [MD_HSS] = HSS2D_compress_nsym(INFO_TREE,params);
        TD_comp = toc; 
        fprintf('HSS-D Compression Time: %e \n',TD_comp)
    
        TD_comp = mean(TD_comp); 
    else
       TD_comp = 0; 
    end
    
    fprintf('\n (I.2) HSS - O(N) \n'); 
    params.Fopts = 'HSS'; params.lev_HSS = lev_HSS; 
    
    tic; 
    [MC_HSS] = HSS2D_compress_nsym(INFO_TREE,params);
    TC_comp = toc; 
    fprintf('HSS-C Compression Time: %e \n',TC_comp)
    
    TC_comp = mean(TC_comp); 
    
else
    if Ntot<Ndense
        fprintf('\n (I.1) Dense - O(Nlog(N)) \n');
        params.Fopts = 'dense';
     
        tic; 
        [MD_HSS] = HSS2D_compress_fsym(INFO_TREE,params);
        TD_comp = toc; 
        fprintf('HSS-D Compression Time: %e \n',TD_comp)
        
    else
       TD_comp = 0; 
    end
    
    TD_comp = mean(TD_comp); 
    
    fprintf('\n (I.2) HSS - O(N) \n'); 
    params.Fopts = 'HSS'; params.lev_HSS = lev_HSS; 
    tic; 
    [MC_HSS] = HSS2D_compress_fsym(INFO_TREE,params);
    TC_comp = toc; 
    fprintf('HSS-C Compression Time: %e \n',TC_comp)
    
    TC_comp = mean(TC_comp); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n (III) Matrix Apply \n'); 
Ntot = size(X,1); 
m = 1; 
apply_params.transinv = TI; 
apply_params.INTERform = params.INTERform; 

% For non-oscillatory problems, we test on m random vectors (gaussian,
% normalized) for small N, and for random columns of the identity. 
if Ntot<10000
   if params.kh == 0 | TI==1
        fprintf('\n Random right hand side \n'); 
        Vec = randn(Ntot,m);    
        Vec = Vec.*(ones(Ntot,1)*(1./sqrt(sum(Vec.*Vec))));
   else
       % The oscillatory problem tested for in this code corresponds to the
        % Lippman Schwinger equation (scattering), for which we build the corresponding
        % right hand side. 
        fprintf('\n Lippman Schwinger, scattered incoming field \n'); 
        Vec = exp(1i*params.kh*X(:,1)); 
        Bump = Bump_function(X(:,1),X(:,2),4); 
        Vec = -params.kh^2*Bump.*Vec; 
   end
   
else
   p = randperm(Ntot); 
   Jp = p(1:m); 
   Id = eye(Ntot,m); 
   Vec(p,:) = Id; 
end

% For the symmetric case, we use the translation invariant apply (apply_TI)
% if TI==1, since it makes full use of vectorization to perform one matrix
% apply per level. 
if params.sym == 0
    if Ntot<Ndense
        apply_params.Fopts = 'dense';
        tic; 
        Y_dense = HSS2D_apply_nsym(INFO_TREE,MD_HSS,Vec,apply_params);
        Td_Apply = toc/m; 
    end
    
    apply_params.Fopts = 'HSS'; apply_params.lev_HSS = lev_HSS; 
    tic; 
    Y_HSS   = HSS2D_apply_nsym(INFO_TREE,MC_HSS,Vec,apply_params);
    TH_Apply = toc/m; 
    
else
    apply_params.Fopts = 'dense';
    if Ntot<Ndense
        tic; 
        if TI == 1 
            Y_dense = HSS2D_apply_TI(INFO_TREE,MD_HSS,Vec,apply_params);
            Td_Apply = toc/m; 
        else
            tic;
            Y_dense = HSS2D_apply_fsym(INFO_TREE,MD_HSS,Vec,apply_params); 
            Td_Apply = toc/m; 
        end
    else
        Td_Apply = 0; 
    end
        
    apply_params.Fopts = 'HSS'; apply_params.lev_HSS = lev_HSS; 
    
    tic; 
    if TI == 1
       Y_HSS   = HSS2D_apply_TI(INFO_TREE,MC_HSS,Vec,apply_params);
    else
       Y_HSS   = HSS2D_apply_fsym(INFO_TREE,MC_HSS,Vec,apply_params); 
    end
    TH_Apply = toc/m;
    
end

if params.dim == 3
   fprintf('\n 3D shape\n');  
   X = shape_3D(X,params,params.type);
   display(params.type)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Relative Error measures 

% If N is small, we compare against the dense apply (forming the entire
% matrix). 
if Ntot<10000       
    A = Kernel_Eval(X,X,params);
    Y_true = A*Vec;
    Errd_Apply = sqrt(max(sum((Y_true-Y_dense).^2)));
    ErrH_Apply = sqrt(max(sum((Y_true-Y_HSS).^2)));
    Err_HSS    = sqrt(max(sum((Y_dense-Y_HSS).^2)));  
    
    fprintf('\n Error dense apply: %e \n',Errd_Apply)
    fprintf('\n Error HSS apply: %e \n',ErrH_Apply)
    
% else, we form the solution to applying the permutation of identity, which
% corresponds to certain columns of our matrix. 
elseif Ntot<Ndense
    Err_HSS    = norm(Y_dense-Y_HSS)/norm(Y_dense);
    
    Y_true = Kernel_Eval(X,X(Jp,:),params); 
    Errd_Apply = sqrt(max(sum((Y_true-Y_dense).^2)));
    ErrH_Apply = sqrt(max(sum((Y_true-Y_HSS).^2)));
    
    fprintf('\n Error dense apply: %e \n',Errd_Apply)
    fprintf('\n Error HSS apply: %e \n',ErrH_Apply)
else
    Errd_Apply = 0;
    ErrH_Apply = 0;
    Err_HSS    = 0;
    
end


Err_DvsC = Err_HSS;
Err_D = Errd_Apply; 
Err_C = ErrH_Apply; 
TDA = Td_Apply; 
TCA = TH_Apply; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Storage measurements 

% Binary tree
info_T = whos('INFO_TREE'); MB_T = info_T.bytes/(2^20); RpDOF_T = info_T.bytes/(8*Ntot); 

% HSS2D Matrices
if Ntot<Ndense
info_AD = whos('MD_HSS'); MB_AD = info_AD.bytes/(2^20) + MB_T; RpDOF_AD = info_AD.bytes/(8*Ntot) + RpDOF_T; 
end
info_AC = whos('MC_HSS'); MB_AC = info_AC.bytes/(2^20) + MB_T; RpDOF_AC = info_AC.bytes/(8*Ntot) + RpDOF_T; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display

% Relative Errors
fprintf('\n Mean Error of HSS vs Dense Matrix Application: %e \n',Err_DvsC);
fprintf('\n Mean O(N^1.5) Apply Time: %e \n',TDA);
fprintf('\n Mean O(N) Apply Time: %e \n',TCA);

% Storage in MB and Reals stored per degree of freedom
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