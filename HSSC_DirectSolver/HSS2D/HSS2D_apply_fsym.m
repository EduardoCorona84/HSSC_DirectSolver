function Y = HSS2D_apply_fsym(INFO_TREE,A_HSS,X,params)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%     FUNCTION CALL:
%         Y = HSS2D_apply_fsym(INFO_TREE,A_HSS,X,params)
%
%     DESCRIPTION:
%
%         This function performs HSS2D matrix-vector multiply AX=Y, where A is symmetric (fsym)
%         and X is a subset of sample points {xi}. The HSS2D form of A is encoded in
%         INFO_TREE, containing the binary tree, index sets, and interpolation operators and
%         A_HSS, containing the self/sibling interaction matrices. We note that INFO_TREE stores
%         L and R and A_HSS stores D in the telescoping factorization.
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
%         are lowrank.
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

trinv = params.transinv; 
num_boxes = length(INFO_TREE.BOX);
phi_up = cell(num_boxes,1);
phi_dn = phi_up;
depth = INFO_TREE.depth;

% Upward Pass ("S2M and M2M")
for lev = depth:-1:1
        if trinv == 1
            k = INFO_TREE.LEV(lev+1).k; 
            J = INFO_TREE.LEV(lev+1).J; 
            T = INFO_TREE.LEV(lev+1).T;
            n = length(J); 
        end
    for i = 1:INFO_TREE.numlev(lev+1)
        Nbox = INFO_TREE.box_numbers(lev+1,i);
        % Index arrays for box Nbox
        I_src = INFO_TREE.BOX(Nbox).I_src;
        if trinv == 0
        k = INFO_TREE.BOX(Nbox).k; 
        J = INFO_TREE.BOX(Nbox).J; 
        T = INFO_TREE.BOX(Nbox).T;
        n = length(J); 
        end
        
        % If the box is a leaf
        if INFO_TREE.BOX(Nbox).child(1) == 0
            q = X(I_src,:); 
        else
            %Children box numbers
            C1 = INFO_TREE.BOX(Nbox).child(1); 
            C2 = INFO_TREE.BOX(Nbox).child(2); 
            q = [phi_up{C1} ; phi_up{C2} ]; 
        end
        
        % Multiply times R
              
        if (k<n)
            if (LR == 0 || isempty(T.V)) %lev>=depth-1)
                phi_up{Nbox} = q(J(1:k),:) + T.U*q(J(k+1:end),:); 
            else
                phi_up{Nbox} = q(J(1:k),:) + T.U*(T.V.'*q(J(k+1:end),:)); 
            end
        else
            phi_up{Nbox} = q(J,:); 
        end
        
        % Upward density
        %phi_up{Nbox} = tmp; 
        %clear tmp
        
    end
end

Y = zeros(size(X)); 


% Downward Pass ("M2L, L2L and L2T")
for lev = 0:depth
        if trinv == 1
            k = INFO_TREE.LEV(lev+1).k; 
            J = INFO_TREE.LEV(lev+1).J; 
            T = INFO_TREE.LEV(lev+1).T;
        end
    for i = 1:INFO_TREE.numlev(lev+1)
        Nbox = INFO_TREE.box_numbers(lev+1,i);
        % Index arrays for box Nbox
        I_src = INFO_TREE.BOX(Nbox).I_src;
        if trinv == 0
            k = INFO_TREE.BOX(Nbox).k; 
            J = INFO_TREE.BOX(Nbox).J; 
            T = INFO_TREE.BOX(Nbox).T;        
        end
        
        % Multiply by D (Sibling Interactions / M2L ) 
        
        if INFO_TREE.BOX(Nbox).child(1) == 0
            if trinv == 0
                Bself = A_HSS(Nbox).Bself; 
            else
                if i==1
                    Bself = A_HSS(lev+1).Bself; 
                end
            end
            
            if lev>lev_HSS
                lam = Bself*X(I_src,:); 
            else
                lam = HSS1D_apply_fsym(Bself,X(I_src,:)); 
            end
        else
            C1 = INFO_TREE.BOX(Nbox).child(1);
            C2 = INFO_TREE.BOX(Nbox).child(2);
            
            if trinv == 0
                k1 = INFO_TREE.BOX(C1).k; 
                B12 = A_HSS(Nbox).B12; 
                B21 = A_HSS(Nbox).B21;
            else
                k1 = INFO_TREE.LEV(INFO_TREE.BOX(C1).levbox + 1).k; 
                B12 = A_HSS(lev+1).B12; 
                B21 = A_HSS(lev+1).B21;
            end
            
            %lam = zeros(k1+k2,1); 
            
            if lev>lev_HSS   
                lam = [ B12*phi_up{C2} ; B21*phi_up{C1} ]; 
            else
                %lam = [ HSS1D_apply_fsym(B12,phi_up{C2}) ; HSS1D_apply_fsym(B21,phi_up{C1}) ]; 
                lam = HSS1D_apply_fsym(B12,[phi_up{C1} ; phi_up{C2}]);
            end
        end
                  
        % For boxes except the top, multiply by L
        if lev == 0
            eta = zeros(size(lam)); 
        else
            eta = zeros(size(lam)); 
            tmp = phi_dn{Nbox}; 
            
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
        %clear tmp
        
        % Downward densities
        phi = eta + lam; 
        clear eta lam; 
        
        % For leaf boxes, phi is Y(I_box)
        if INFO_TREE.BOX(Nbox).child(1) == 0
            Y(I_src,:) = phi; 
        else
        %Children box numbers 
            phi_dn{C1} = phi(1:k1,:); 
            phi_dn{C2} = phi(k1+1:end,:); 
        end
        
    end
end