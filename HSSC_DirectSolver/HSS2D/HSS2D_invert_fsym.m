function A_HSS = HSS2D_invert_fsym(INFO_TREE,params)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%     FUNCTION CALL:
%         A_HSS = HSS2D_invert_fsym(INFO_TREE,params)
%
%     DESCRIPTION:
%         General HSS inversion is explained in section 2.4 of the paper (Algorithm 1), where
%         the matrices D,L,R in the factorization of A^-1 are given in terms of Schur complement
%         matrices E and F^-1. This function computes E and F^-1 for each box with the
%         assumption that A is symmetric (fsym).
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
%         INFO_TREE   <struct>    (output of HSS2D_build_tree_fsym) Binary tree information and
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
%             A_HSS(i).E          (see HSS2D_build_E_fsym).
%





Fopts_def = params.Fopts;
Fopts = 'dense';
Eopts.mat = 'dense'; 
if strcmp(params.INTERform,'dense')
    EoptsLR_def = 0; 
else
    EoptsLR_def = 1; 
end

TI = params.transinv;

if ~strcmp(Fopts_def,'dense') 
    lev_HSS = params.lev_HSS; 
else
    lev_HSS = -1;
end

fprintf('\n levHSS = %d ',lev_HSS); 
rho = zeros(INFO_TREE.depth+1,1);
current_lev = INFO_TREE.depth;  

% Upward Pass
fprintf('\n lev = '); 
for lev = current_lev:-1:0  
    fprintf(' %d ,',lev);   
    
    % If kernel is Translation Invariant (TI), only one (k,J,T) per level
    % is stored in INFO_TREE.LEV(lev+1)
    if TI == 1
        k = INFO_TREE.LEV(lev+1).k;
        J = INFO_TREE.LEV(lev+1).J;
        T = INFO_TREE.LEV(lev+1).T;
        numlev = 1; %only one box per level
    else
        numlev = INFO_TREE.numlev(lev+1); 
    end
        
    for i = 1:numlev
        Nbox = INFO_TREE.box_numbers(lev+1,i);
        
        %------------------------------------------------------------------
        % Routine to get blocks of FInv from E
        
        % If kernel is Non Translation Invariant (NTI), one (k,J,T) per box
        % is stored in INFO_TREE.BOX(Nbox)
        if TI == 0
            k = INFO_TREE.BOX(Nbox).k;
            J = INFO_TREE.BOX(Nbox).J;
            T = INFO_TREE.BOX(Nbox).T;
        end
           
        % If the box is a leaf
        if INFO_TREE.BOX(Nbox).child(1) == 0
            X = params.X_source; 
            % Index arrays for box Nbox
            I_src = INFO_TREE.BOX(Nbox).I_src;
            
            if TI == 0
                A_HSS(Nbox).FI = inv(Kernel_Eval(X(I_src,:),X(I_src,:),params)); 
            else
                A_HSS(lev+1).FI = inv(Kernel_Eval(X(I_src,:),X(I_src,:),params)); 
            end
        else
            %Children box numbers
            C1 = INFO_TREE.BOX(Nbox).child(1); C2 = INFO_TREE.BOX(Nbox).child(2); 
            % Skeleton indices for c1 and c2
            I1_sk = INFO_TREE.BOX(C1).I_sk; 
            I2_sk = INFO_TREE.BOX(C2).I_sk; 
            
            % Build FInv: 
            
            %N = length(I_src);   
            if (lev>lev_HSS)
                Fopts = 'dense'; 
            elseif (lev == lev_HSS &&  (~strcmp(Fopts_def,'dense') && ~strcmp(Fopts_def,'block')))     
                Fopts = 'd2HSS';
            elseif (lev == lev_HSS && strcmp(Fopts_def,'block'))
                Fopts = 'd2block';
            else
                Fopts = Fopts_def; 
            end
            
            if ~strcmp(Fopts_def,'dense')
                if TI == 0
                    params.JE1 = A_HSS(C1).JE; params.JE2 = A_HSS(C2).JE;  
                else
                    params.JE1 = A_HSS(lev+2).JE; params.JE2 = A_HSS(lev+2).JE; 
                end
            end
            
            params.rho = rho; 
            
            if TI == 0
                %Build FInv
                [FI] = HSS2D_build_Finv_fsym( A_HSS(C1).E,    A_HSS(C2).E,    lev+1, X, J, k, I1_sk,I2_sk,params,Fopts); 
                A_HSS(Nbox).FI = FI; 
            else
                %Build FInv
                [FI] = HSS2D_build_Finv_fsym( A_HSS(lev+2).E, A_HSS(lev+2).E, lev+1, X, J, k, I1_sk,I2_sk,params,Fopts); 
                A_HSS(lev+1).FI = FI; 
            end
            
            clear FI; 
        end
        %------------------------------------------------------------------
        
        %Routine to compute E from FInv
        if lev>0 || (strcmp(params.type,'torus') && params.m<1)      

            if (strcmp(Fopts_def,'block') && INFO_TREE.BOX(Nbox).child(1) == 0)
                Eopts.mat = 'dense'; 
            elseif strcmp(Fopts,'d2HSS') 
                Eopts.mat = 'HSS';
            elseif strcmp(Fopts,'d2block')
                Eopts.mat = 'block'; 
            else
                Eopts.mat = Fopts;
                %if (~strcmp(Fopts,'dense') && ~strcmp(Fopts,'block'))  
                %    Eopts.mat = 'HSS'; 
                %end
            end
            
            if lev < INFO_TREE.depth-1
                Eopts.LR = EoptsLR_def; 
            else
                Eopts.LR = 0;
            end
        
            if TI == 0
                [A_HSS(Nbox).E,A_HSS(Nbox).JE] = HSS2D_build_E_fsym(A_HSS(Nbox).FI, J, T, k, params.acc, Eopts);
            else
                [A_HSS(lev+1).E,A_HSS(lev+1).JE] = HSS2D_build_E_fsym(A_HSS(lev+1).FI, J, T, k, params.acc, Eopts);
            end
            
            if (strcmp(Fopts_def,'HSS') || strcmp(Fopts_def,'HSStest')) && lev<=lev_HSS
                if TI == 0
                    rho(lev+1) = max(rho(lev+1),size(A_HSS(Nbox).FI.rs2sk.U,2)); 
                else
                    rho(lev+1) = size(A_HSS(lev+1).FI.rs2sk.U,2);
                end
            end
        end
        %------------------------------------------------------------------
        
    end
end

end
