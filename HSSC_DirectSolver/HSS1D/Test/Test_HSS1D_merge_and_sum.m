function Test_HSS1D_merge_and_sum(TREE,params)
%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%   Test code testing the interaction of: 
%   
%   - HSS1D merge of two diagonal sub-blocks.
%   - Off-diagonal Compression of an HSS matrix given a binary tree
%   - Sum of the two HSS matrices resulting from the previous two
%   
%   This mimics the HSS2D_buildFInv routine. 
%   
%   1) We take the I_sk arrays from a uniform binary tree, sort them and compress 
%   A_sk = K[I_sk,I_sk] for 2 pairs of sibling boxes. We use HSS1D_split to
%   separate skeleton and interface blocks, and then merge diagonal blocks of each 
%   pair. 
%   
%   2) We then compress the off-diagonal blocks on an HSS structure, and use
%   HSS1D sum to add this to the diagonal blocks. 
%   
%   Finally, we test these by applying the corresponding matrices to N random
%   vectors and measuring timings and error. 
%

global BOX 
rand( 'seed',0);

% 2D TREE indices and parameters
X = params.X_source;
h = params.h; lay = params.layers; acc = params.acc; 

% Test 1: Sibling boxes on even level (2)
Nbox1 = 4;
Nbox2 = 5;
Nbox3 = 6; 
Nbox4 = 7; 

tic; toc; 

lev = TREE.BOX(7).levbox;
I1_sk = TREE.BOX(Nbox1).I_sk; X1_sk = X(I1_sk,:); 
I2_sk = TREE.BOX(Nbox2).I_sk; X2_sk = X(I2_sk,:);  
I3_sk = TREE.BOX(Nbox3).I_sk; X3_sk = X(I3_sk,:); 
I4_sk = TREE.BOX(Nbox4).I_sk; X4_sk = X(I4_sk,:);

cB1 = TREE.BOX(Nbox1).cent;
cB2 = TREE.BOX(Nbox2).cent;
cB3 = TREE.BOX(Nbox3).cent;
cB4 = TREE.BOX(Nbox4).cent;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (I) Find interface sets between Box 1 and 2, and between Box 3 and 4. 
fprintf('\n(I) Find Interface Sets 12 and 34 \n') 
I1 = LOCAL_Find_Interface(X1_sk,cB1,lev,h,lay,'west');
I2 = LOCAL_Find_Interface(X2_sk,cB2,lev,h,lay,'south');
I3 = LOCAL_Find_Interface(X3_sk,cB3,lev,h,lay,'north');
I4 = LOCAL_Find_Interface(X4_sk,cB4,lev,h,lay,'east');

n1 = length(I1_sk); n2 = length(I2_sk); 
n3 = length(I3_sk); n4 = length(I4_sk);

% Binary trees to perform compression 
T1 = LOCAL_get_tree(n1,56);
T2 = LOCAL_get_tree(n2,56);
T3 = LOCAL_get_tree(n3,56);
T4 = LOCAL_get_tree(n4,56);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (II) Compress HSS matrices
fprintf('\n(II) Compress HSS Matrices given trees \n') 
if params.sym == 0
    tic;
    %{
    E1 = HSS1D_compress_nsym_giventree(X1_sk,T1,1e-10,'Lap',params,1);
    E2 = HSS1D_compress_nsym_giventree(X2_sk,T2,1e-10,'Lap',params,1);
    E3 = HSS1D_compress_nsym_giventree(X3_sk,T3,1e-10,'Lap',params,1);
    E4 = HSS1D_compress_nsym_giventree(X4_sk,T4,1e-10,'Lap',params,1);
    %}
    E1 = HSS1D_compress_nsym('green',X1_sk,T1,1e-10,params);
    E2 = HSS1D_compress_nsym('green',X2_sk,T2,1e-10,params);
    E3 = HSS1D_compress_nsym('green',X3_sk,T3,1e-10,params);
    E4 = HSS1D_compress_nsym('green',X4_sk,T4,1e-10,params);
    T_comp_sk = toc/4; 
else
    tic; 
    %{
    E1 = HSS1D_compress_fsym_giventree(X1_sk,T1,1e-10,'Lap',params,1);
    E2 = HSS1D_compress_fsym_giventree(X2_sk,T2,1e-10,'Lap',params,1);
    E3 = HSS1D_compress_fsym_giventree(X3_sk,T3,1e-10,'Lap',params,1);
    E4 = HSS1D_compress_fsym_giventree(X4_sk,T4,1e-10,'Lap',params,1);
    %}
    E1 = HSS1D_compress_fsym('green',X1_sk,T1,1e-10,params,1);
    E2 = HSS1D_compress_fsym('green',X2_sk,T2,1e-10,params,1);
    E3 = HSS1D_compress_fsym('green',X3_sk,T3,1e-10,params,1);
    E4 = HSS1D_compress_fsym('green',X4_sk,T4,1e-10,params,1);
    T_comp_sk = toc/4; 
end

fprintf('(II) Compression Time %2.5f \n',T_comp_sk) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (III) SPLIT
fprintf('\n(III) Split Into Boundary / Interface Sets \n') 
if params.sym == 0
    tic;
    [E1_sk,E1_rs] = HSS1D_split_nsym(E1,I1,I1,acc); 
    [E2_sk,E2_rs] = HSS1D_split_nsym(E2,I2,I2,acc); 
    [E3_sk,~]     = HSS1D_split_nsym(E3,I3,I3,acc); 
    [E4_sk,~]     = HSS1D_split_nsym(E4,I4,I4,acc); 
    T_split = toc/4; 
else
    tic;
    [E1_sk,E1_rs] = HSS1D_split_fsym(E1,I1,acc); 
    [E2_sk,E2_rs] = HSS1D_split_fsym(E2,I2,acc); 
    [E3_sk,~]     = HSS1D_split_fsym(E3,I3,acc); 
    [E4_sk,~]     = HSS1D_split_fsym(E4,I4,acc); 
    T_split = toc/4; 
end

fprintf('(III) Split Time %2.5f \n',T_split) 

XI_sk  = [X1_sk(I1,:) ; X2_sk(I2,:)]; 
XI_rs  = [X1_sk(~I1,:) ; X2_sk(~I2,:)];
XII_sk = [X3_sk(I3,:) ; X4_sk(I4,:)]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (IV) MERGE 
fprintf('\n(IV) Merge Diagonal Blocks 12 and 34 \n') 
if params.sym == 0
    tic;
    fprintf('\nMerge 12\n') 
    [E12_sk]  = HSS1D_merge_nsym(E1_sk,E2_sk,XI_sk,lev,acc,'concatenate'); 
    fprintf('\nMerge 34\n')
    [E34_sk]  = HSS1D_merge_nsym(E3_sk,E4_sk,XII_sk,lev,acc,'concatenate'); 
    T_merge = toc/2; 
    [E12_rs,X12rs,Jrs]  = HSS1D_merge_nsym(E1_rs,E2_rs,XI_rs,lev,acc,'align');
else
    tic;
    fprintf('\nMerge 12\n') 
    [E12_sk]  = HSS1D_merge_fsym(E1_sk,E2_sk,XI_sk,lev,acc,'concatenate');  
    fprintf('\nMerge 34\n')
    [E34_sk] = HSS1D_merge_fsym(E3_sk,E4_sk,XII_sk,lev,acc,'concatenate'); 
    T_merge = toc/2; 
    [E12_rs,X12rs,Jrs]  = HSS1D_merge_fsym(E1_rs,E2_rs,XI_rs,lev,acc,'align'); 
end

fprintf('(IV) Merge Time %2.5f \n',T_merge) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k1 = sum(I1); k2 = sum(I2); k3 = sum(I3); k4 = sum(I4);

% Dense matrices
KI_skdg  = [Kernel_Eval(X1_sk(I1,:),X1_sk(I1,:),params) zeros(k1,k2) ; zeros(k2,k1) Kernel_Eval(X2_sk(I2,:),X2_sk(I2,:),params)]; 
KI_skoff = [zeros(k1,k1) Kernel_Eval(X1_sk(I1,:),X2_sk(I2,:),params) ; Kernel_Eval(X2_sk(I2,:),X1_sk(I1,:),params) zeros(k2,k2)]; 
KI_sk = Kernel_Eval(XI_sk,XI_sk,params); 

r1 = sum(~I1); r2 = sum(~I2); 
KI_rsdg  = [Kernel_Eval(X1_sk(~I1,:),X1_sk(~I1,:),params) zeros(r1,r2) ; zeros(r2,r1) Kernel_Eval(X2_sk(~I2,:),X2_sk(~I2,:),params)]; 
KI_rsoff = [zeros(r1,r1) Kernel_Eval(X1_sk(~I1,:),X2_sk(~I2,:),params) ; Kernel_Eval(X2_sk(~I2,:),X1_sk(~I1,:),params) zeros(r2,r2)]; 
KI_rs = Kernel_Eval(XI_rs,XI_rs,params); 

KI_rsdg = KI_rsdg(Jrs,Jrs); KI_rsoff = KI_rsoff(Jrs,Jrs); KI_rs = KI_rs(Jrs,Jrs); 

KII_skdg = [Kernel_Eval(X3_sk(I3,:),X3_sk(I3,:),params) zeros(k3,k4) ; zeros(k4,k3) Kernel_Eval(X4_sk(I4,:),X4_sk(I4,:),params)];
KII_skoff = [zeros(k3,k3) Kernel_Eval(X3_sk(I3,:),X4_sk(I4,:),params) ; Kernel_Eval(X4_sk(I4,:),X3_sk(I3,:),params) zeros(k4,k4)]; 
KII_sk = Kernel_Eval(XII_sk,XII_sk,params); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Off-diagonal compression and sum
fprintf('\n(V) OffDiagonal Compression and Sum\n') 
if params.sym == 0
    %E12_off  = HSS1D_compress_nsym_offd_giventree(XI_sk,E12_sk,1e-12,'lap',params,1); 
    E12_off  = HSS1D_compress_nsym('green',XI_sk,E12_sk,1e-12,params,'offd'); 
    E12_sum  = HSS1D_sum_nsym(E12_sk,E12_off,1e-12); 
    
    Brs = [true(1,sum(~I1)) false(1,sum(~I2))]; 
    Brs = Brs(Jrs); 
    E12_rsoff  = HSS1D_compress_nsym_offd_leaf(X12rs,E12_rs,Brs,1e-12,params); 
    E12_rssum  = HSS1D_sum_nsym(E12_rs,E12_rsoff,acc);

    %E34_off = HSS1D_compress_nsym_offd_giventree(XII_sk,E34_sk,1e-12,'lap',params,1); 
    E34_off  = HSS1D_compress_nsym('green',XII_sk,E34_sk,1e-12,params,'offd'); 
    E34_sum  = HSS1D_sum_nsym(E34_sk,E34_off,1e-12); 
else
    E12_off  = HSS1D_compress_fsym('green',XI_sk,E12_sk,1e-12,params,'offd'); 
    E12_sum  = HSS1D_sum_fsym(E12_sk,E12_off,acc); 
    
    Brs = [true(1,sum(~I1)) false(1,sum(~I2))]; 
    Brs = Brs(Jrs); 
    E12_rsoff  = HSS1D_compress_fsym_offd_leaf(X12rs,E12_rs,Brs,1e-12,params); 
    E12_rssum  = HSS1D_sum_fsym(E12_rs,E12_rsoff,acc); 

    E34_off = HSS1D_compress_fsym('green',XII_sk,E34_sk,1e-12,params,'offd'); 
    E34_sum  = HSS1D_sum_fsym(E34_sk,E34_off,acc); 
end

M = 10; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n(VI) Apply Experiments\n') 
Vec12 = randn(k1+k2,M); 
Vecrs12 = randn(r1+r2,M); 
Vec34 = randn(k3+k4,M); 
    
if params.sym == 0
    U12_merge    = HSS1D_apply_nsym(E12_sk,Vec12,1); 
    U12_mergesum = HSS1D_apply_nsym(E12_sum,Vec12,1); 
    U12HSS_offd  = HSS1D_apply_nsym(E12_off,Vec12,1); 
    
    U12rs_merge    = HSS1D_apply_nsym(E12_rs,Vecrs12); 
    U12rs_mergesum = HSS1D_apply_nsym(E12_rssum,Vecrs12); 
    U12rsHSS_offd  = HSS1D_apply_nsym(E12_rsoff,Vecrs12);
else
    U12_merge    = HSS1D_apply_fsym(E12_sk,Vec12); 
    U12_mergesum = HSS1D_apply_fsym(E12_sum,Vec12); 
    U12HSS_offd  = HSS1D_apply_fsym(E12_off,Vec12); 
        
    U12rs_merge    = HSS1D_apply_fsym(E12_rs,Vecrs12); 
    U12rs_mergesum = HSS1D_apply_fsym(E12_rssum,Vecrs12); 
    U12rsHSS_offd  = HSS1D_apply_fsym(E12_rsoff,Vecrs12); 
end
    
U12      = KI_sk*Vec12; 
U12_dg   = KI_skdg*Vec12; 
U12_offd = KI_skoff*Vec12; 
    
U12rs      = KI_rs*Vecrs12; 
U12_rsdg   = KI_rsdg*Vecrs12; 
U12_rsoffd = KI_rsoff*Vecrs12; 
    
if params.sym == 0 
    U34_merge    = HSS1D_apply_nsym(E34_sk,Vec34,1); 
    U34_mergesum = HSS1D_apply_nsym(E34_sum,Vec34,1); 
    U34HSS_offd  = HSS1D_apply_nsym(E34_off,Vec34,1); 
else
    U34_merge    = HSS1D_apply_fsym(E34_sk,Vec34); 
    U34_mergesum = HSS1D_apply_fsym(E34_sum,Vec34); 
    U34HSS_offd  = HSS1D_apply_fsym(E34_off,Vec34); 
end
    
U34     = KII_sk*Vec34; 
U34_dg   = KII_skdg*Vec34; 
U34_offd = KII_skoff*Vec34; 
   
Err_offd     = max(norm(U12_offd - U12HSS_offd) / norm(U12_offd),norm(U34_offd - U34HSS_offd) / norm(U34_offd)); 
Err_dg       = max(norm(U12_dg - U12_merge) / norm(U12_dg),norm(U34_dg - U34_merge) / norm(U34_dg)); 
Err_mergesum = max(norm(U12 - U12_mergesum) / norm(U12),norm(U34 - U34_mergesum) / norm(U34)); 
    
Err_rsoffd     = norm(U12_rsoffd - U12rsHSS_offd) / norm(U12_rsoffd); 
Err_rsdg       = norm(U12_rsdg - U12rs_merge) / norm(U12_rsdg); 
Err_rsmergesum = norm(U12rs - U12rs_mergesum) / norm(U12rs); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n Merge Error: \n') 
mean(Err_dg)

fprintf('\n Off Diagonal Compression Error: \n') 
mean(Err_offd)

fprintf('\n Sum Error: \n') 
mean(Err_mergesum)

fprintf('\n RS Merge Error: \n') 
mean(Err_rsdg)

fprintf('\n RS Off Diagonal Compression Error: \n') 
mean(Err_rsoffd)

fprintf('\n RS Sum Error: \n') 
mean(Err_rsmergesum)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NODES = LOCAL_get_tree(ntot,nbox_max)
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

nboxes = ibox_last;
NODES  = cell(50,nboxes);
for ibox = 1:nboxes
  for j = 1:7
    NODES{j,ibox} = TREE{j,ibox};
  end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I_bd = LOCAL_Find_Interface(X,cB,lev,h,lay,type)

% Center points X
CB = repmat(cB,size(X,1),1);
XC = X - CB;

% Dimensions
X1max = max(XC(:,1)); 
X2max = max(XC(:,2));
Lsize = lay*h - h/2;

k = size(X,1); 

% Separate in North, South, East and West
if mod(lev,2) == 1
    IN = (XC(:,2)>=(X2max - Lsize));
    IS = (XC(:,2)<=-(X2max - Lsize));
    IE = (abs(XC(:,2))<(X2max - Lsize) & (XC(:,1)>=(X1max - Lsize)));
    IW = (abs(XC(:,2))<(X2max - Lsize) & (XC(:,1)<=-(X1max - Lsize)));
else
    IN = (XC(:,2)>=(X2max - Lsize)) & (abs(XC(:,1))<(X1max - Lsize));
    IS = (XC(:,2)<=-(X2max - Lsize)) & (abs(XC(:,1))<(X1max - Lsize));
    IE = (XC(:,1)>=(X1max - Lsize));
    IW = (XC(:,1)<=-(X1max - Lsize));
end

I_bd = true(1,k); 

if mod(lev,2) == 0
    if strcmp(type,'north') || strcmp(type,'west')
        I_bd(IS) = false; 
    elseif strcmp(type,'east') || strcmp(type,'south')
        I_bd(IN) = false; 
    end
else
    if strcmp(type,'north') || strcmp(type,'west')
        I_bd(IE) = false; 
    elseif strcmp(type,'south') || strcmp(type,'east')
        I_bd(IW) = false; 
    end
end

end