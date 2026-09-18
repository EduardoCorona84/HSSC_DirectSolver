function K_HSS = HSS1D_compressLR_fsym(U,V,TREE,acc)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        K_HSS = HSS1D_compressLR_fsym(U,V,TREE,acc)
%
%    DESCRIPTION:
%        This function compresses the low-rank fully-symmetric matrix UV' into HSS form (K_HSS).
%        It is used for compressing the interpolation operators in HSS2D.
%
%    INFILE VARIABLES:
%        radius_rel (double) specifies the relative radius of the proxy circles. default = 1.75.
%        nproxy (int) specifies the number of samples to be used on the proxy circles. default = 50.
%
%    INPUT:
%        U is the left low rank (n x q) matrix
%        V is the right low rank (n x q) matrix
%        TREE is the {7 x nboxes} cell array that gives the underlying binary tree structure of
%            K_HSS.
%        acc (double) is the desired accuracy.
%
%    OUTPUT:
%        K_HSS is the {50 x nboxes} cell array of the HSS matrix UV'.
%


global BOX

% Compute the tree strcture.
ntot   = size(U,1);

if length(TREE) == 1
    nbox_max = TREE; 
    K_HSS  = LOCAL_get_tree(ntot,nbox_max);
else 
    nboxes = size(TREE,2); 
    K_HSS  = cell(50,nboxes);
    K_HSS(BOX.DATA,:) = TREE(BOX.DATA,:); 
    clear TREE;
end

nboxes = size(K_HSS,2);

% Compress all K_HSS, going from smaller to larger.
for ibox = nboxes:(-1):2

  % Construct the index vectors for ibox.
  if ( (K_HSS{BOX.C1,ibox}<=0) && (K_HSS{BOX.C2,ibox}<=0) ) % ibox has no sons.
      indskel = K_HSS{BOX.END1,ibox} - 1 + (1:K_HSS{BOX.END2,ibox});
  else % ibox has two sons
    ison1       = K_HSS{BOX.C1,ibox};
    ison2       = K_HSS{BOX.C2,ibox};
    indskel = [K_HSS{BOX.I_SKUP,ison1},K_HSS{BOX.I_SKUP,ison2}];
  end
  
  A21 = V(indskel,:).';  
        
  % Compute the skeletons.
  [T,J] = ID(A21, acc);
  k = size(T,1);
  
  % Record the outgoing skeletons:
  K_HSS{BOX.I_SKUP,ibox} = indskel(J(1:k));
  K_HSS{BOX.T_UP,ibox} = T;
  K_HSS{BOX.J_UP,ibox} = J;
  
  % Record nskel and kskel
  K_HSS{BOX.NSKEL,ibox} = length(indskel);
  K_HSS{BOX.KSKEL,ibox} = k;
end

% Construct the matrices representing self-interactions on the leaves.
for ibox = nboxes:(-1):2
  if ( (K_HSS{BOX.C1,ibox}<=0) && (K_HSS{BOX.C2,ibox}<=0) )
    ind = K_HSS{BOX.END1,ibox} - 1 + (1:K_HSS{BOX.END2,ibox});
    
    K_HSS{BOX.M_SELF,ibox} = U(ind,:)*(V(ind,:).'); 
    
  end
end

% Construct the matrices for sibling interactions.
for ibox = nboxes:(-1):1
  if ( (K_HSS{BOX.C1,ibox}>0) && (K_HSS{BOX.C2,ibox}>0) )
    ison1 = K_HSS{BOX.C1,ibox};
    ison2 = K_HSS{BOX.C2,ibox};
    
    K_HSS{BOX.M_SIB,ison1} = U(K_HSS{BOX.I_SKUP,ison1},:)*(V(K_HSS{BOX.I_SKUP,ison2},:).');
    K_HSS{BOX.M_SIB,ison2} = K_HSS{BOX.M_SIB,ison1}.'; 
    
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function K_HSS = LOCAL_get_tree(ntot,nbox_max)
global BOX 

TREE = cell(7,ntot);

% Construct the top node.
TREE{BOX.LEVEL,1} = 0;
TREE{BOX.PARENT,1} = NaN;
TREE{BOX.C1,1} = -1;
TREE{BOX.C2,1} = -1;
TREE{BOX.END1,1} = 1;
TREE{BOX.END2,1} = ntot;
ibox_last = 0;
ncreated  = 1;
ilevel    = 0;

% Create smaller K_HSS via hierarchical subdivision.
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

nboxes = ibox_last;
K_HSS  = cell(50,nboxes);
for ibox = 1:nboxes
  for j = 1:7
    K_HSS{j,ibox} = TREE{j,ibox};
  end
end

return