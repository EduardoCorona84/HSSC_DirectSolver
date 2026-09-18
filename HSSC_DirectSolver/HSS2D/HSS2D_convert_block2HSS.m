function AH = HSS2D_convert_block2HSS(Ab,TREE,params)
%{
 This function is mainly for testing purposes: 

Inputs: 

Ab - produced by HSS2D_Inverse_Compression with params.Fopts = 'block' 
(modified inversion algorithm, produces E and blocks of F in dense form). 
TREE - 2D binary tree + skeleton information
params

Output
AH - Brute-force compression of top blocks in HSS1D form (analogous to
output of HSS2D_Inverse_Compression with params.Fopts = 'HSS')

One can change 'brute' for 'green' on the compression routines. In that
case, kernel evaluation parameters should be passed as an argument. 
%}

AH = Ab; 
acc = params.acc; 

%if params.sym == 1
%    TREE = HSS2D_build_tree_fsym(params); 
%else
%    TREE = HSS2D_build_tree_nsym(params); 
%end

lev_HSS = levHSS(TREE,params.n_cut,params.transinv); 
display(lev_HSS)

for lev = lev_HSS:-1:0
    
    if params.transinv == 0
        num = TREE.numlev(lev+1); 
    else
        num = 1; 
    end
    
    for i = 1:num 
        if params.transinv == 0
           ibox = TREE.box_numbers(lev+1,i); 
        else
           ibox = lev+1;  
        end
        
        % Matrix sizes
        k = size(Ab(ibox).FI.SrsInv,2); 
        r = size(Ab(ibox).FI.rs_Inv,2);
        
        % HSS1D binary trees
        Trs = LOCAL_tree(ceil(r/2),floor(r/2)); 
        Tsk = LOCAL_tree(ceil(k/2),floor(k/2)); 
        
        if params.sym == 1
            % F_rs^{-1} and S_rs^{-1} HSS1D compression (F_r2s and Fs2r are the same). 
            AH(ibox).FI.rs_Inv  = HSS1D_compress_fsym('brute',Ab(ibox).FI.rs_Inv,Trs,acc); 
            AH(ibox).FI.SrsInv = HSS1D_compress_fsym('brute', Ab(ibox).FI.SrsInv,Tsk,acc); 
            if ibox>1
                % E HSS1D compression
                AH(ibox).E = HSS1D_compress_fsym('brute',Ab(ibox).E,Tsk,acc); 
            end
        
        else
            % F_rs^{-1} and S_rs^{-1} HSS1D compression (F_r2s and Fs2r are the same). 
            AH(ibox).FI.rs_Inv  = HSS1D_compress_nsym('brute',Ab(ibox).FI.rs_Inv,Trs,acc); 
            AH(ibox).FI.SrsInv = HSS1D_compress_nsym('brute', Ab(ibox).FI.SrsInv,Tsk,acc); 
            if ibox>1
                % E HSS1D compression
                AH(ibox).E = HSS1D_compress_nsym('brute',Ab(ibox).E,Tsk,acc); 
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverse Apply to compare block with HSS-block Inverses
fprintf('\n Inverse Apply \n'); 
X = params.X_source;
Ntot = size(X,1); 
m = 10; 

% random rhs
Vec = randn(Ntot,m); 
Vec = Vec.*(ones(Ntot,1)*(1./sqrt(sum(Vec.*Vec))));

params.lev_HSS = lev_HSS; params.Fopts = 'HSS';
if params.sym == 1
    if params.transinv == 1
        Y_HSS   = HSS2D_applyinv_TI(TREE,AH,Vec,params);
        params.Fopts = 'block'; 
        Y_block = HSS2D_applyinv_TI(TREE,Ab,Vec,params);
    else
        Y_HSS   = HSS2D_applyinv_fsym(TREE,AH,Vec,params);
        params.Fopts = 'block'; 
        Y_block = HSS2D_applyinv_fsym(TREE,Ab,Vec,params);
    end
else
    Y_HSS   = HSS2D_applyinv_nsym(TREE,AH,Vec,params);
    params.Fopts = 'block'; 
    Y_block = HSS2D_applyinv_nsym(TREE,Ab,Vec,params);
end
    


if Ntot<2000
    A = Kernel_Eval(X,X,params);
    Y_true = A\Vec;
    ErrH_Apply = sqrt(max(sum((Y_true-Y_HSS).^2)));
    Err_HSS    = sqrt(max(sum((Y_block-Y_HSS).^2))); 
    Errb_Apply  = sqrt(max(sum((Y_true - Y_block).^2))); 
else
    params.Fopts = 'dense'; 
    if sym == 1
        [M_HSS] = HSS2D_compress_fsym(TREE,params);
        VH      = HSS2D_apply_fsym(TREE,M_HSS,Y_HSS,params);
        VB      = HSS2D_apply_fsym(TREE,M_HSS,Y_block,params);
    else
        [M_HSS] = HSS2D_compress_nsym(TREE,params);
        VH      = HSS2D_apply_nsym(TREE,M_HSS,Y_HSS,params);
        VB      = HSS2D_apply_nsym(TREE,M_HSS,Y_block,params);
    end
    
    
    ErrH_Apply = sqrt(max(sum((Vec-VH).^2)));
    Err_HSS    = sqrt(max(sum((Y_block-Y_HSS).^2))); 
    Errb_Apply = sqrt(max(sum((Vec-VB).^2))); 
end

mEH = Err_HSS; 
mEHA = ErrH_Apply; 

fprintf('\n Mean Error of HSS vs Dense Block Inverse Application: %e \n',mEH);
fprintf('\n Mean Error of HSS Inverse Application: %e \n',mEHA);
fprintf('\n Mean Error of Dense Block Inverse Application: %e \n',Errb_Apply);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TREE = LOCAL_tree(n1,n2)
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

% Split in half
%nbox_max = ceil(min(n1,n2)/2);
nbox_max = 56; 
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