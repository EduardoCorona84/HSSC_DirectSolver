function Y = HSS2D_applyinv_TI(INFO_TREE,A_HSS,X,params)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%     FUNCTION CALL:
%         Y = HSS2D_applyinv_TI(INFO_TREE,A_HSS,X,params)
%
%     DESCRIPTION:
%         For transation-invariant kernels (TI), this function performs HSS2D inverse matrix-
%         vector multiply Y = A^{-1}X, where X is a subset of sample points. The HSS2D form of A
%         is encoded in INFO_TREE, containing the binary tree, index sets, and interpolation
%         operators and A_HSS, containing the Schur complement matrices E and F^-1 of each box.
%         We note that INFO_TREE stores the telescoping factorization components L,R of A and
%         A_HSS stores L~,D~,R~ of A^-1.
%
%         This is an implementation of Algorithm 6 in the paper. It follows the general HSS
%         matvec algorithm (section 2.3), using L~,D~,R~ of the telescoping factorization of
%         A^-1. The only difference is that L~,D~,R~ must first be computed from E,F^-1,L,R at
%         each box,
%
%                 R~ = ERF^1
%                 D~ = F^-1(I - LR~)
%                 L~ = F^-1LE,
%
%          where these operations respect that E,L,R are lowrank and F^-1 is partitioned into
%          2x2 compressed blocks. Additional vectorization speedup results from the translation
%          invariance, where the same set of matrices is applied for every box in a given level.
%
%
%     INPUT:
%         INFO_TREE   <struct>    (output of HSS2D_build_tree_fsym) HSS binary tree and matrix
%                                     skeletons / interpolation matrices of A.
%         A_HSS       <struct>    (output of HSS2D_invert_fsym) Schur complement
%                                     matrices E,F^-1 used to compute the L~,D~,R~ components of
%                                     the telescoping factorization of A^-1.
%         X           <Nx2>       Column vector of charges on sample points.
%         params      <struct>    (output of HSS_tree_parameters) Kernel parameters.
%
%
%     OUTPUT:
%         Y           <Nx1>       Result of the HSS2D inverse-matrix apply Y = A^{-1}X.
%




Fopts_def = params.Fopts;
if strcmp(Fopts_def,'HSS') || strcmp(Fopts_def,'block')
    lev_HSS = params.lev_HSS; 
else
    lev_HSS = -1; 
end

if strcmp(params.INTERform,'dense')
    LR = 0;
else
    LR = 1; 
end

m = size(X,2); 

Y = zeros(size(X)); 
%nboxes = length(INFO_TREE.BOX); 
depth = INFO_TREE.depth; 
phi_up = cell(depth+1,1); q = phi_up; phi_dn = phi_up; 

% Upward Pass ("S2M and M2M")
for lev = depth:-1:1
    % Determine Fopts and Eopts based on level
    if lev>lev_HSS
        Fopts = 'dense';
        Eopts = 'dense'; 
    else
        Fopts = Fopts_def; 
        if strcmp(Fopts,'block')
            Eopts = 'dense'; 
        else
            Eopts = Fopts_def; 
        end
    end
    
    k = INFO_TREE.LEV(lev+1).k;
    J = INFO_TREE.LEV(lev+1).J;
    T = INFO_TREE.LEV(lev+1).T;
    n = length(J); 
    FI = A_HSS(lev+1).FI; 
    E = A_HSS(lev+1).E; 
    JE = A_HSS(lev+1).JE; 
    
    nblev = INFO_TREE.numlev(lev+1);
    q{lev+1} = zeros(n,nblev*m); 
    tmp = zeros(k,nblev*m); phi_up{lev+1} = tmp; 
        
    for i = 1:nblev
        Nbox = INFO_TREE.box_numbers(lev+1,i);
        % Index arrays for box Nbox
        I_src = INFO_TREE.BOX(Nbox).I_src;
        
        % If the box is a leaf
        if lev == depth
            q{lev+1}(:,m*(i-1)+1:m*i) = X(I_src,:); 
        else
            %Children box numbers
            Nt = (2^(lev+1)-1);
            C1 = INFO_TREE.BOX(Nbox).child(1); i1 = C1-Nt;
            C2 = INFO_TREE.BOX(Nbox).child(2); i2 = C2-Nt;
            q{lev+1}(:,m*(i-1)+1:m*i) = [phi_up{lev+2}(:,m*(i1-1)+1:m*i1) ; phi_up{lev+2}(:,m*(i2-1)+1:m*i2) ]; 
        end
    end

    psi = HSS2D_applyFinv_fsym(FI,q{lev+1},k,Fopts);
        
    if (k<n)
        if (LR == 0 || isempty(T.V)) %lev>=depth-1)
            tmp = psi(J(1:k),:) + T.U*psi(J(k+1:end),:); 
        else
            tmp = psi(J(1:k),:) + T.U*(T.V.'*psi(J(k+1:end),:)); 
        end
    else
        tmp = psi(J,:); 
    end
        
    % Upward density (Apply E) 
        
    if strcmp(Eopts,'dense')
        phi_up{lev+1}(JE,:) = E*tmp(JE,:); 
    elseif strcmp(Eopts,'HSS')
        phi_up{lev+1}(JE,:) = HSS1D_apply_fsym(E,tmp(JE,:));  
    end
        
    clear k J T FI E JE tmp psi; 
end

% Downward Pass ("M2L, L2L and L2T")
for lev = 0:depth
    % Determine Fopts and Eopts based on level
    if lev>lev_HSS
        Fopts = 'dense';
        Eopts = 'dense'; 
    else
        Fopts = Fopts_def; 
        if strcmp(Fopts,'block')
            Eopts = 'dense'; 
        else
            Eopts = Fopts_def; 
        end
    end
    
    k = INFO_TREE.LEV(lev+1).k; 
    J = INFO_TREE.LEV(lev+1).J; 
    T = INFO_TREE.LEV(lev+1).T;
    n = length(J); 
    FI = A_HSS(lev+1).FI; 
    E = A_HSS(lev+1).E; 
    JE = A_HSS(lev+1).JE; 
    nblev = INFO_TREE.numlev(lev+1);
    
    if lev == 0
        lam = [phi_up{2}(:,1:m) ; phi_up{2}(:,m+1:end)];
        eta = zeros(size(lam)); 
    else
        tmp = zeros(n,m*nblev);
        if (k < n)
            if (LR == 0 || isempty(T.V)) %lev>=depth-1)  
                tmp(J,:) = [phi_up{lev+1} ; T.U.'*phi_up{lev+1}]; 
            else
                tmp(J,:) = [phi_up{lev+1} ; (T.V)*(T.U.'*phi_up{lev+1}) ];
            end
        else
            tmp(J,:) = phi_up{lev+1}; 
        end
            
        lam = q{lev+1} - tmp;
        
        % For boxes except the top, multiply by Lhat = FInv*L*E
        eta = zeros(size(lam)); tmp = zeros(k,m*nblev); 
          
        % Apply E
        if strcmp(Eopts,'dense')
            tmp(JE,:) = E*phi_dn{lev+1}(JE,:); 
        elseif strcmp(Eopts,'HSS')
            tmp(JE,:) = HSS1D_apply_fsym(E,phi_dn{lev+1}(JE,:));  
        end 
            
        if (k<n)
            if (LR == 0 || lev>=depth-1)
                eta(J,:) = [tmp ; T.U.'*tmp ]; 
            else
                eta(J,:) = [tmp ; (T.V)*(T.U.'*tmp) ]; 
            end
        else
            eta(J,:) = tmp; 
        end
    end
    clear tmp; 
    
    % Downward densities
    phi = HSS2D_applyFinv_fsym(FI,eta + lam,k,Fopts);
    clear eta lam k J T FI E JE;   
        
    if lev<depth
        % Number of children boxes
        Nc = 2^(lev+1);        
        kc = INFO_TREE.LEV(lev+2).k; 
        
        M  = repmat((1:m),2^lev,1)'; M = M(:)'; 
        I1 = repmat((1:2:Nc),m,1); I1 = I1(:)'; 
        I2 = repmat((2:2:Nc),m,1); I2 = I2(:)';
        
        I1 = m*(I1-1)+M; I2 = m*(I2-1)+M;
        
        phi_dn{lev+2}(:,I1) = phi(1:kc,:); 
        phi_dn{lev+2}(:,I2) = phi(kc+1:end,:); 
    else
        for i = 1:nblev
        Nbox = INFO_TREE.box_numbers(depth+1,i);
        % Index arrays for box Nbox
        I_src = INFO_TREE.BOX(Nbox).I_src; 
        Y(I_src,:) = phi(:,m*(i-1)+1:m*i); 
        end
    end
    
end