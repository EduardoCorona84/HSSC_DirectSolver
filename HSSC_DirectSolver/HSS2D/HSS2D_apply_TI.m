function Y = HSS2D_apply_TI(INFO_TREE,A_HSS,X,params)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%     FUNCTION CALL:
%         Y = HSS2D_apply_TI(INFO_TREE,A_HSS,X,params)
%
%     DESCRIPTION:
%         For translation invariant kernels (TI), this function performs HSS2D matrix-vector
%         multiply AX=Y, where X is a subset of sample points {xi}. The HSS2D form of A is
%         encoded in INFO_TREE, containing the binary tree, index sets, and interpolation
%         operators and A_HSS, containing the self/sibling interaction matrices. We note that
%         INFO_TREE stores L and R and A_HSS stores D in the telescoping factorization.
%
%         HSS matrix-vector multiplication is outlined in section 2.3 of the paper, and it
%         follows directly from the recursive version of the telescoping factorization. The
%         product is computed by keeping track of two sets of vectors, x^l and u^l (in the
%         paper, phi is used in place of x). The vectors x^l represent the charges assigned to
%         the skeleton points at each level in the tree, and they are computed in an upward pass
%
%                     x^d = X
%                     x^l = R^{l+1} x^{l+1}.
%
%         The vectors u^l are the resulting potentials, computed in a downward pass
%
%                     u^0 = 0
%                     u^l = D^l x^l + L^l u^{l-1}
%                     u^d = Y.
%
%         O(N) complexity is achieved in HSS2D because the D matrices are HSS1D and L,R matrices
%         are lowrank. Additional vectorization speedup results from the translation invariance,
%         where the same set of matrices is applied for every box in a given level.
%
%     INPUT:
%         INFO_TREE   <struct>    (output of HSS2D_build_tree_fsym) Binary tree information and
%                                     matrix skeletons / interpolation matrices.
%         A_HSS       <struct>    (output of HSS2D_compress_fsym) Self and sibling
%                                     interaction matrices that make up D in the TF.
%         X           <Nx1 float> Column vector of charges on sample points.
%         params      <struct>    (output of HSS_tree_parameters) Kernel parameters.
%
%
%     OUTPUT:
%         Y   <Nx1 float>         Result of the HSS2D matrix apply Y = AX.
%




% TODO: Have to initialize phi_up, phi_dn and Y. 

opts = params.Fopts;
if strcmp(opts,'HSS') || strcmp(opts,'block')
    lev_HSS = params.lev_HSS; 
else
    lev_HSS = -1; 
end

if strcmp(params.INTERform,'dense')
    LR = 0;
else
    LR = 1; 
end

num_boxes = length(INFO_TREE.BOX);
phi_up = cell(num_boxes,1);
phi_dn = phi_up;
depth = INFO_TREE.depth;
m = size(X,2); 

% Upward Pass ("S2M and M2M")
for lev = depth:-1:0
    k = INFO_TREE.LEV(lev+1).k; 
    J = INFO_TREE.LEV(lev+1).J; 
    T = INFO_TREE.LEV(lev+1).T;
    n = length(J); 
    
    numlev = INFO_TREE.numlev(lev+1);
    q{lev+1} = zeros(n,numlev*m); 
        
    for i = 1:numlev
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
        
    % Multiply times R          
    if (k<n)
        if (LR == 0 || isempty(T.V))% lev>=depth-1)
            phi_up{lev+1} = q{lev+1}(J(1:k),:) + T.U*q{lev+1}(J(k+1:end),:); 
        else
            phi_up{lev+1} = q{lev+1}(J(1:k),:) + T.U*(T.V.'*q{lev+1}(J(k+1:end),:)); 
        end
    else
        phi_up{lev+1} = q{lev+1}(J,:); 
    end
end

Y = zeros(size(X)); 


% Downward Pass ("M2L, L2L and L2T")
for lev = 0:depth
    k = INFO_TREE.LEV(lev+1).k; 
    J = INFO_TREE.LEV(lev+1).J; 
    T = INFO_TREE.LEV(lev+1).T;
    
    if lev == depth 
       if depth>lev_HSS
           lam = A_HSS(lev+1).Bself*q{depth+1};
       else
           lam = HSS1D_apply_fsym(A_HSS(lev+1).Bself,q{depth+1});
       end
    else
        %lam = zeros(n,m*numlev);
        B12 = A_HSS(lev+1).B12; 
        B21 = A_HSS(lev+1).B21;
        
        if lev>lev_HSS
            %lam = [ B12*phi_up{C2} ; B21*phi_up{C1} ]; 
            B = [zeros(size(B12,1),size(B12,1)) B12 ; B21 zeros(size(B21,1),size(B21,1))];
            lam = B*q{lev+1}; 
        else
            lam = HSS1D_apply_fsym(B12,q{lev+1});           
        end
    end
        
    if lev == 0
        eta = zeros(size(lam)); 
    else
        eta = zeros(size(lam)); 
        tmp = phi_dn{lev+1}; 
            
        if (k<n)
            if (LR == 0 || isempty(T.V)) %lev>=depth-1)
                if ~isempty(T.U)
                    eta(J,:) = [tmp ; T.U.'*tmp ]; 
                else
                    eta(J,:) = tmp; 
                end
            else
                eta(J,:) = [tmp ; (T.V)*(T.U.'*tmp) ]; 
            end
        else
            eta(J,:) = tmp; 
        end
    end
        
    % Downward densities
    phi = eta + lam; 
    clear eta lam;   
    
    numlev = INFO_TREE.numlev(lev+1);
    
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
        for i = 1:numlev
        Nbox = INFO_TREE.box_numbers(depth+1,i);
        % Index arrays for box Nbox
        I_src = INFO_TREE.BOX(Nbox).I_src; 
        Y(I_src,:) = phi(:,m*(i-1)+1:m*i); 
        end
    end
end