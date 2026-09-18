function [A_HSS] = HSS2D_compress_fsym(INFO_TREE,params)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%     FUNCTION CALL:
%         [A_HSS] = HSS2D_compress_fsym(INFO_TREE,params)
%
%     DESCRIPTION:
%         Given tree information and parameters, this function returns an array of structs
%         containing the matrices that comprise D in the telescoping factorization. That is, the
%         self interaction matrices for leaf boxes and sibling interaction matrices for non-leaf
%         boxes. It is assumed that the underlying matrix A is symmetric (fsym).
%
%     INPUT:
%         INFO_TREE   <struct>    (output of HSS2D_build_tree_fsym) Binary tree.
%         params      <struct>    (output of HSS_tree_parameters) With added parameters below:
%             .Fopts              {'dense', 'HSS'}
%                                 - 'dense' Matrices at all boxes in the tree are stored
%                                     densely.
%                                 - 'HSS' Matrices at boxes above level .lev_HSS in the tree are
%                                     compressed as HSS1D. Required for the O(N) solver.
%             .lev_HSS            (output of levHSS) Matrices are stored densely (rather than
%                                     HSS1D) on and below this level of the tree. Automatically
%                                     set to -1 if Fopts = 'dense'.
%
%             .T_source and .X_source (IGNORE FOR PLANAR SURFACES) are changed only if the
%                 domain is on a 3D surface. .T_source must be set to the source points used in
%                 the BUILD TREE stage (the old .X_source). .X_Source must be redefined as the
%                 parameterization onto the 3D surface. (Note this is still in the testing phase)
%
%                 .T_source   <Nx2 float>     = params.X_source
%                 .X_source   <Nx3 float>     = [(x(u,v), y(u,v), z(u,v)) for (u,v) in .T_source]
%
%     OUTPUT:
%         A_HSS           <mx1 struct>    D matrices for each box in the tree (or each level for
%                                             TI kernels)
%         A_HSS(i)                        Info for box i (or level i).
%             .Bself      <nxn float>     Dense self interactions matrix for leaves (between
%                                             source points within the box).
%             .B12/.B21   <struct>        Sibling interaction matrix for non-leaves (between
%                                             children skeleton sets X1,X2). It is stored either
%                                             densely or as HSS1D depending on params.Fopts.
%
%                                         B12 = K[X1,X2] , B21 = K[X2,X1] = B12'
%




Source = params.X_source; 
D = params.dim;
if D==3
    display(D)   
   T_src = params.T_source;  
end
opts = params.Fopts;
if strcmp(opts,'HSS') || strcmp(opts,'block')
    lev_HSS = params.lev_HSS; 
else
    lev_HSS = -1; 
end

% Kernel Evaluation parameters
par = params; par.X_source = []; 

TI = params.transinv; 
nmax = 100; 

% Downward Pass 
for lev = 0:INFO_TREE.depth
    if TI == 0
            numlev = INFO_TREE.numlev(lev+1); 
    else
            numlev = 1; 
    end
    
    for i = 1:numlev
        Nbox = INFO_TREE.box_numbers(lev+1,i);
        % Index arrays for box Nbox
        I_src = INFO_TREE.BOX(Nbox).I_src;   
        
        % Multiply by D (Sibling Interactions / M2L ) 
        
        if INFO_TREE.BOX(Nbox).child(1) == 0
            % Self interaction matrix 
            if lev>lev_HSS
                Bself = Kernel_Eval(Source(I_src,:),Source(I_src,:),par); 
            else
                if D == 3
                par.T = T_src(I_src,:); 
                end
                Bself  = HSS1D_compress_fsym('green',Source(I_src,:),nmax,params.acc,par); 
            end
            
            if TI == 0
                A_HSS(Nbox).Bself = Bself;
            else
                A_HSS(lev+1).Bself = Bself;
            end
        else
            C1 = INFO_TREE.BOX(Nbox).child(1);
            I1 = INFO_TREE.BOX(C1).I_sk;
            
            C2 = INFO_TREE.BOX(Nbox).child(2);
            I2 = INFO_TREE.BOX(C2).I_sk;
            
            %lam = zeros(k1+k2,1); 
            X1 = Source(I1,:); X2 = Source(I2,:);
            
            % Sibling Interaction Matrices Evaluation
            if lev>lev_HSS
                B12 = Kernel_Eval(X1,X2,params); 
                B21 = B12.';
            else
                X12 = [X1 ; X2];
                T12 = LOCAL_HSS1D_tree(size(X1,1),size(X2,1),nmax); 
                
                if D == 3
                par.T = [T_src(I1,:) ; T_src(I2,:)];      
                end
                
                B12  = HSS1D_compress_fsym('green',X12,T12,params.acc,par,'offd'); 
                B21 = [];
            end
            
            if TI == 0
                A_HSS(Nbox).B12 = B12; 
                A_HSS(Nbox).B21 = B21; 
            else
                A_HSS(lev+1).B12 = B12; 
                A_HSS(lev+1).B21 = B21; 
            end
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TREE = LOCAL_HSS1D_tree(n1,n2,nbox_max)
global BOX 

ntot = n1+n2; 
TREE = cell(50,3);

% Construct the top node.
TREE{BOX.LEVEL,1} = 0;
TREE{BOX.PARENT,1} = NaN;
TREE{BOX.C1,1} = 2;
TREE{BOX.C2,1} = 3;
TREE{BOX.END1,1} = 1;
TREE{BOX.END2,1} = ntot; 

% Nodes 2 and 3
TREE{BOX.LEVEL,2} = 1; TREE{BOX.LEVEL,3} = 1; 
TREE{BOX.PARENT,2} = 1; TREE{BOX.PARENT,3} = 1;  
TREE{BOX.END1,2} = 1; TREE{BOX.END2,2} = n1; 
TREE{BOX.END1,3} = n1+1; TREE{BOX.END2,3} = n2; 

TREE{BOX.C1,2} = -1; TREE{BOX.C2,2} = -1;
TREE{BOX.C1,3} = -1; TREE{BOX.C2,3} = -1;

% parameters to keep subdividing the tree, if necessary
ibox_last = 1;
ncreated  = 2;
ilevel    = 1;

% Create smaller nodes via hierarchical subdivision.
% We sweep one level at a time.
while (ncreated > 0)
  ibox_first = ibox_last + 1;
  ibox_last  = ibox_last + ncreated;
  ncreated   = 0;
  for ibox = ibox_first:ibox_last
    nbox = TREE{BOX.END2,ibox};
    if (nbox > nbox_max)
      nhalf             = ceil(nbox/2);
      ibox_son1         = ibox_last + ncreated + 1;
      ibox_son2         = ibox_last + ncreated + 2;
      TREE{BOX.C1,ibox}      = ibox_son1;
      TREE{BOX.C2,ibox}      = ibox_son2;
      TREE{BOX.LEVEL,ibox_son1} = ilevel+1;
      TREE{BOX.LEVEL,ibox_son2} = ilevel+1;
      TREE{BOX.PARENT,ibox_son1} = ibox;
      TREE{BOX.PARENT,ibox_son2} = ibox;
      TREE{BOX.C1,ibox_son1} = -1;
      TREE{BOX.C1,ibox_son2} = -1;
      TREE{BOX.C2,ibox_son1} = -1;
      TREE{BOX.C2,ibox_son2} = -1;
      TREE{BOX.END1,ibox_son1} = TREE{BOX.END1,ibox};
      TREE{BOX.END1,ibox_son2} = TREE{BOX.END1,ibox} + nhalf;
      TREE{BOX.END2,ibox_son1} = nhalf;
      TREE{BOX.END2,ibox_son2} = TREE{BOX.END2,ibox} - nhalf;
      ncreated          = ncreated + 2;
    end
  end
  ilevel = ilevel + 1;
end

end