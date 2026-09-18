function A_HSS = HSS2D_invert_nsym(INFO_TREE,params)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%     FUNCTION CALL:
%         A_HSS = HSS2D_invert_nsym(INFO_TREE,params)
%
%     DESCRIPTION:
%         General HSS inversion is explained in section 2.4 of the paper (Algorithm 1), where
%         the matrices D,L,R in the factorization of A^-1 are given in terms of Schur complement
%         matrices E and F^-1. This function computes E and F^-1 for each box without making any
%         assumptons on the symmetry of A (nsym).
%
%         Computation of E and F^-1 at each box requires inverting an nxn matrix which takes
%         O(N^3/2) time for HSS2D (params.Fopts = 'dense'). To achieve linear complexity we
%         store/compute interpolation operators as low-rank, E as HSS1D, and blocks of F^-1 as
%         HSS1D or low-rank (params.Fopts = 'HSS'). This is outlined in Algorithm 2, and
%         detailed in Algorithms 3,4,5.
%
%         For translation-invariant kernels, one set of matrices per level in the tree need to
%         be computed, rather than per box.
%
%     INPUT:
%         INFO_TREE   <struct>    (output of HSS2D_build_tree_nsym) Binary tree information and
%                                     L,R matrices in the telescoping factorization.
%         params      <struct>    (output of HSS_tree_parameters) With added parameters below.
%             .Fopts              {'dense', 'block', 'HSS'}
%                                 - 'dense' F^-1 and E are each stored and computed as a single
%                                     dense matrix : O(N^3/2).
%                                 - 'block' E and blocks of F are stored densely : O(N^3/2).
%                                 - 'HSS' blocks of E and F are compressed : O(N).
%             .lev_HSS            (output of levHSS) Matrices are stored densely on and below
%                                     this level of the tree. Automatically set to -1 if Fopts =
%                                     'dense'.
%
%     OUTPUT:
%         A_HSS   <mx1 struct>    A_HSS(i) contains the matrices E and F^-1 used to represent
%                                 the inverse of A.
%             A_HSS(i).FI.*       (see HSS2D_build_Finv_fsym).
%             A_HSS(i).E          (see HSS2D_build_E_nsym).
%





Fopts_def = params.Fopts;
Fopts = 'dense';
Eopts.mat = 'dense'; 
if strcmp(params.INTERform,'dense')
    EoptsLR_def = 0; 
else
    EoptsLR_def = 1; 
end

if ~strcmp(Fopts_def,'dense') 
    lev_HSS = params.lev_HSS; 
else
    lev_HSS = -1;
end

display(lev_HSS)

rho = zeros(INFO_TREE.depth+1,1); 

% Upward Pass
for lev = INFO_TREE.depth:-1:0        
    for i = 1:INFO_TREE.numlev(lev+1)
        Nbox = INFO_TREE.box_numbers(lev+1,i);
        
        %------------------------------------------------------------------
        % Routine to get blocks of FInv from E
        
        k = INFO_TREE.BOX(Nbox).k;
        J = INFO_TREE.BOX(Nbox).J;
        T_up = INFO_TREE.BOX(Nbox).T_up;
        T_dn = INFO_TREE.BOX(Nbox).T_dn;
        
        % If the box is a leaf
        if INFO_TREE.BOX(Nbox).child(1) == 0
            X = params.X_source; 
            % Index arrays for box Nbox
            I_src = INFO_TREE.BOX(Nbox).I_src;
            
            A_HSS(Nbox).FI = inv(Kernel_Eval(X(I_src,:),X(I_src,:),params)); 
        else
            %Children box numbers
            C1 = INFO_TREE.BOX(Nbox).child(1); 
            C2 = INFO_TREE.BOX(Nbox).child(2); 
            % Skeleton indices for c1 and c2
            I1_sk = INFO_TREE.BOX(C1).I_sk; 
            I2_sk = INFO_TREE.BOX(C2).I_sk; 
            
            % Build FInv: 
            
            if (lev>lev_HSS)
                Fopts = 'dense'; 
            elseif (lev == lev_HSS && (strcmp(Fopts_def,'HSS') || strcmp(Fopts_def,'HSStest')))
                Fopts = 'd2HSS';
            elseif (lev == lev_HSS && strcmp(Fopts_def,'block'))
                Fopts = 'd2block';
            else
                Fopts = Fopts_def; 
            end
            
            if ~strcmp(Fopts_def,'dense')
                params.JE1 = A_HSS(C1).JE; params.JE2 = A_HSS(C2).JE;  
            end
            
            params.rho = rho; 
            
            %Build FInv
            [FI] = HSS2D_build_Finv_nsym( A_HSS(C1).E,    A_HSS(C2).E,    lev+1, X, J, k, I1_sk,I2_sk,params,Fopts); 
            A_HSS(Nbox).FI = FI; 
            
            clear FI; 
        end
        %------------------------------------------------------------------
        
        %Routine to compute E from FInv
        if lev>0

            if (strcmp(Fopts_def,'block') && INFO_TREE.BOX(Nbox).child(1) == 0)
                Eopts.mat = 'dense'; 
            elseif strcmp(Fopts,'d2HSS') || strcmp(Fopts,'d2HSSb')
                Eopts.mat = 'HSS';
            elseif strcmp(Fopts,'d2block')
                Eopts.mat = 'block'; 
            else
                Eopts.mat = Fopts;
                if strcmp(Fopts,'HSS') || strcmp(Fopts,'HSStest')
                    Eopts.mat = 'HSS'; 
                end
            end
            
            if lev < INFO_TREE.depth-1
                Eopts.LR = EoptsLR_def; 
            else
                Eopts.LR = 0;
            end
        
            % Build E
            [A_HSS(Nbox).E,A_HSS(Nbox).JE] = HSS2D_build_E_nsym(A_HSS(Nbox).FI, J, T_up, T_dn, k, params.acc, Eopts);
            
            if (strcmp(Fopts_def,'HSS') || strcmp(Fopts_def,'HSStest')) && lev<=lev_HSS
                rho(lev+1) = max(rho(lev+1),size(A_HSS(Nbox).FI.rs2sk.U,2)); 
            end
        end
        %------------------------------------------------------------------
        
    end
end

end
