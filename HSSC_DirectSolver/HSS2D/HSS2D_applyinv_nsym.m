function Y = HSS2D_applyinv_nsym(INFO_TREE,A_HSS,X,params)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%     FUNCTION CALL:
%         Y = HSS2D_applyinv_nsym(INFO_TREE,A_HSS,X,params)
%
%     DESCRIPTION:
%         Without making any assumptions on the symmetry of the underlying kernel, this function
%         performs HSS2D inverse matrix-vector multiply Y = A^{-1}X, where X is a subset of
%         sample points. The HSS2D form of A is encoded in INFO_TREE, containing the binary
%         tree, index sets, and interpolation operators and A_HSS, containing the Schur
%         complement matrices E and F^-1 of each box. We note that INFO_TREE stores the
%         telescoping factorization components L,R of A and A_HSS stores L~,D~,R~ of A^-1.
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
%          2x2 compressed blocks.
%
%      INPUT:
%         INFO_TREE   <struct>    (output of HSS2D_build_tree_nsym) HSS binary tree and matrix
%                                     skeletons / interpolation matrices of A.
%         A_HSS       <struct>    (output of HSS2D_invert_nsym) Schur complement
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
if strcmp(Fopts_def,'HSStest')
    Fopts_def = 'HSS'; 
end

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
nboxes = length(INFO_TREE.BOX); 
phi_up = cell(nboxes,1); phi_dn = phi_up; 
transinv = params.transinv; 
depth = INFO_TREE.depth; 

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
    
    nblev = INFO_TREE.numlev(lev+1);
        
    for i = 1:nblev
        Nbox = INFO_TREE.box_numbers(lev+1,i);
        % Index arrays for box Nbox
        I_src = INFO_TREE.BOX(Nbox).I_src;
        
        k = INFO_TREE.BOX(Nbox).k;
        J = INFO_TREE.BOX(Nbox).J; 
        T = INFO_TREE.BOX(Nbox).T_up;
        n = length(J); 
        FI = A_HSS(Nbox).FI; 
        E = A_HSS(Nbox).E; 
        JE = A_HSS(Nbox).JE;
        
        % If the box is a leaf
        if INFO_TREE.BOX(Nbox).child(1) == 0
            q = X(I_src,:); 
        else
            %Children box numbers
            C1 = INFO_TREE.BOX(Nbox).child(1); 
            C2 = INFO_TREE.BOX(Nbox).child(2); 
            q = [phi_up{C1} ; phi_up{C2} ]; 
        end

        psi = HSS2D_applyFinv_nsym(FI,q,k,Fopts);
        
        if (k<n)
            if (LR == 0 || lev>=depth-1)
                tmp = psi(J(1:k),:) + T.U*psi(J(k+1:end),:); 
            else
                tmp = psi(J(1:k),:) + T.U*(T.V.'*psi(J(k+1:end),:)); 
            end
        else
            tmp = psi(J,:); 
        end
        
        % Upward density (Apply E) 
        
        if strcmp(Eopts,'dense')
            phi_up{Nbox}(JE,:) = E*tmp(JE,:); 
        elseif strcmp(Eopts,'HSS')
            phi_up{Nbox}(JE,:) = HSS1D_apply_nsym(E,tmp(JE,:));  
        end
        clear tmp;
    end
    clear k J T FI E JE
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
    
    nblev = INFO_TREE.numlev(lev+1);
    
    for i = 1:nblev
        Nbox = INFO_TREE.box_numbers(lev+1,i);
        % Index arrays for box Nbox
        I_src = INFO_TREE.BOX(Nbox).I_src;
        k = INFO_TREE.BOX(Nbox).k; 
        J = INFO_TREE.BOX(Nbox).J;
        T = INFO_TREE.BOX(Nbox).T_dn;
        n = length(J); 
        FI = A_HSS(Nbox).FI; 
        E = A_HSS(Nbox).E; 
        JE = A_HSS(Nbox).JE;
        
        phi_up_Nbox = phi_up{Nbox};
        
        % Multiply by Di = FInv(I - LRhat) (Sibling Interactions / M2L )  
        tmp = zeros(n,m);
        if lev == 0
            lam = [phi_up{2} ; phi_up{3} ];
        else
            if (k < n)
                if (LR == 0 || lev>=depth-1)
                    tmp(J,:) = [phi_up_Nbox ; T.U.'*phi_up_Nbox ]; 
                else
                    tmp(J,:) = [phi_up_Nbox ; (T.V)*(T.U.'*phi_up_Nbox) ];
                end
            else
                tmp(J,:) = phi_up_Nbox; 
            end
        
            if INFO_TREE.BOX(Nbox).child(1) == 0
                q = X(I_src,:); 
            else
                %Children box numbers
                C1 = INFO_TREE.BOX(Nbox).child(1); 
                C2 = INFO_TREE.BOX(Nbox).child(2); 
                q = [phi_up{C1} ; phi_up{C2} ]; 
            end
            
            lam = q - tmp;
        end
        clear tmp q
        
        % For boxes except the top, multiply by Lhat = FInv*L*E
        if lev == 0
            eta = zeros(size(lam)); 
        else
            eta = zeros(size(lam)); 
            
            % Apply E
            if strcmp(Eopts,'dense')
                tmp(JE,:) = E*phi_dn{Nbox}(JE,:); 
            elseif strcmp(Eopts,'HSS')
                tmp(JE,:) = HSS1D_apply_nsym(E,phi_dn{Nbox}(JE,:));  
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
        clear tmp
        
        % Downward densities
        phi = HSS2D_applyFinv_nsym(FI,eta + lam,k,Fopts);
        clear eta lam; 
        
        % For leaf boxes, phi is Y(I_box)
        if INFO_TREE.BOX(Nbox).child(1) == 0
            Y(I_src,:) = phi; 
        else
        %Children box numbers
            C1 = INFO_TREE.BOX(Nbox).child(1); 
            C2 = INFO_TREE.BOX(Nbox).child(2);
            if transinv == 0
                k1 = INFO_TREE.BOX(C1).k; 
            else
                k1 = INFO_TREE.LEV(lev+2).k; 
            end
            phi_dn{C1} = phi(1:k1,:); 
            phi_dn{C2} = phi(k1+1:end,:); 
        end
        
    end
    clear k J T FI E JE; 
end