function [X_sort,J,I_bd,TREE,card] = sort_skeleton(X,cB,lev,h,lay,type)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%     FUNCTION CALL:
%         [X_sort,J,I_bd,TREE,card] = sort_skeleton(X,cB,lev,h,lay,type)
%
%     DESCRIPTION:
%         This function sorts points on a box's boundary layers, clockwise around the box
%         center. It is generally called on the source points of a box and it is used so that
%         the proxy-skeleton interaction matrix can be compressed as HSS1D.
%
%     INPUT:
%         X       <nx2 float>     Source points that we want sorted.
%         cB      <1x2 float>     Center of the box.
%         lev     <int>           Level in the tree.
%         h       <float>         Width/length of one grid cell.
%         lay     <int>           Number of boundary layers in the skeleton set.
%         type    <string>        {'north', east', 'south', west'} Edge that comes first.
%
%     OUTPUT:
%         X_sort  <kx2 float>     Sorted points: clockwise edges around the box.
%         J       <kx1 int>       Permutation vector such that X_sort = X(J,:)
%         I_bd    <kx1 bool>      Array that is true for boundary points.
%         TREE    <struct>        Restricted tree (HSS1D) for boundary points.
%         card    <1x4 int>       Sizes for north,east,south,west sides of X.
%




% Center points X
CB = repmat(cB,size(X,1),1);
XC = X - CB;

% Dimensions
X1max = max(XC(:,1)); 
X2max = max(XC(:,2));
Lsize = (lay-1)*h + h/2;


k = size(X,1); 
J = 1:k; 

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

% Sort North and South (according to X1), East and West (according to X2)
[~,IN_sort] = sort(XC(IN,1));
[~,IS_sort] = sort(-XC(IS,1));
[~,IE_sort] = sort(-XC(IE,2));
[~,IW_sort] = sort(XC(IW,2));
nN = length(IN_sort); nE = length(IE_sort); nS = length(IS_sort); nW = length(IW_sort); 

% Indices by cardinal point
JN = J(IN); JE = J(IE); JS = J(IS); JW = J(IW); 

% Sort Indices
% Vertical Interfase - [NE E S W NW]

if mod(lev,2) == 0
    if strcmp(type,'north')
        J = [JN(IN_sort) JE(IE_sort) JS(IS_sort) JW(IW_sort)];
        I_bd = [true(size(IN_sort')) true(size(IE_sort')) false(size(IS_sort')) true(size(IW_sort'))];
        
        TREE = LOCAL_get_tree(nN, nE, ceil(nW/2), floor(nW/2)); 
        card = [nN nE nS nW];
    elseif strcmp(type,'east')
        J = [JE(IE_sort) JS(IS_sort) JW(IW_sort) JN(IN_sort)];
        I_bd = [true(size(IE_sort')) true(size(IS_sort')) true(size(IW_sort')) false(size(IN_sort'))];
        
        TREE = LOCAL_get_tree(nE, nS, ceil(nW/2), floor(nW/2));
        card = [nE nS nW nN];
    elseif strcmp(type,'south')
        J = [JS(IS_sort) JW(IW_sort) JN(IN_sort) JE(IE_sort)];
        I_bd = [true(size(IS_sort')) true(size(IW_sort')) false(size(IN_sort')) true(size(IE_sort')) ];
        
        TREE = LOCAL_get_tree(nS, nW, ceil(nE/2), floor(nE/2));
        card = [nS nW nN nE];
    elseif strcmp(type,'west')
        J = [JW(IW_sort) JN(IN_sort) JE(IE_sort) JS(IS_sort)];
        I_bd = [true(size(IW_sort')) true(size(IN_sort')) true(size(IE_sort')) false(size(IS_sort'))];
        
        TREE = LOCAL_get_tree(nW, nN, ceil(nE/2), floor(nE/2));
        card = [nW nN nE nS];
    end
else
    if strcmp(type,'north')
        J = [JN(IN_sort) JE(IE_sort) JS(IS_sort) JW(IW_sort)];
        I_bd = [true(size(IN_sort')) false(size(IE_sort')) true(size(IS_sort')) true(size(IW_sort'))];
        
        TREE = LOCAL_get_tree(ceil(nN/2), floor(nN/2), nS, nW);
        card = [nN nE nS nW];
    elseif strcmp(type,'south')
        J = [JS(IS_sort) JW(IW_sort) JN(IN_sort) JE(IE_sort)];
        I_bd = [true(size(IS_sort')) false(size(IW_sort')) true(size(IN_sort')) true(size(IE_sort')) ];
        
        TREE = LOCAL_get_tree(ceil(nS/2), floor(nS/2), nN, nE);
        card = [nS nW nN nE];
    elseif strcmp(type,'east')
        J = [JE(IE_sort) JS(IS_sort) JW(IW_sort) JN(IN_sort)];
        I_bd = [true(size(IE_sort')) true(size(IS_sort')) false(size(IW_sort')) true(size(IN_sort'))];
        
        TREE = LOCAL_get_tree(nE, nS, ceil(nN/2), floor(nN/2));    
        card = [nE nS nW nN];
    elseif strcmp(type,'west')
        J = [JW(IW_sort) JN(IN_sort) JE(IE_sort) JS(IS_sort)];
        I_bd = [true(size(IW_sort')) true(size(IN_sort')) false(size(IE_sort')) true(size(IS_sort'))];
        
        TREE = LOCAL_get_tree(nW, nN, ceil(nS/2), floor(nS/2));  
        card = [nW nN nE nS];
    end
end

% Sort X
X_sort = X(J,:); 

return

function NODES = LOCAL_get_tree(n1,n2,n3,n4)
global BOX 

ntot = n1+n2+n3+n4; 
TREE = cell(7,7);

% Construct the top node.
TREE{BOX.LEVEL,1} = 0;
TREE{BOX.PARENT,1} = NaN;
TREE{BOX.C1,1} = 2;
TREE{BOX.C2,1} = 3;
TREE{BOX.END1,1} = 1;
TREE{BOX.END2,1} = ntot;

% Construct first two levels
%n1 = length(I1); n2 = length(I2); n3 = length(I3); n4 = length(I4); 

% Nodes 2 and 3
TREE{BOX.LEVEL,2} = 1; TREE{BOX.LEVEL,3} = 1; 
TREE{BOX.PARENT,2} = 1; TREE{BOX.PARENT,3} = 1; 
TREE{BOX.C1,2} = 4; TREE{BOX.C2,2} = 5; TREE{BOX.C1,3} = 6; TREE{BOX.C2,3} = 7; 
TREE{BOX.END1,2} = 1; TREE{BOX.END2,2} = n1+n2; 
TREE{BOX.END1,3} = TREE{BOX.END2,2}+1; TREE{BOX.END2,3} = n3+n4; 

% Nodes 4,5,6,7
TREE{BOX.LEVEL,4} = 2; TREE{BOX.LEVEL,5} = 2; TREE{BOX.LEVEL,6} = 2; TREE{BOX.LEVEL,7} = 2; 
TREE{BOX.PARENT,4} = 2; TREE{BOX.PARENT,5} = 2; TREE{BOX.PARENT,6} = 3; TREE{BOX.PARENT,7} = 3; 
for i=4:7
TREE{BOX.C1,i} = -1; TREE{BOX.C2,i} = -1;
end

TREE{BOX.END2,4} = n1; TREE{BOX.END2,5} = n2;   TREE{BOX.END2,6} = n3;      TREE{BOX.END2,7} = n4; 

fullb = zeros(1,7); 
fullb(1) = 1; 
nz = 0; 
for ibox = 1:3
    ison1 = TREE{BOX.C1,ibox}; ison2 = TREE{BOX.C2,ibox};  
   if TREE{BOX.END2,ison1} > 0
      fullb(ison1) = 1;
      TREE{BOX.C1,ibox} = TREE{BOX.C1,ibox} - nz; 
   else
       TREE{BOX.C1,ibox} = -1; 
       nz = nz+1; 
   end
   
   if TREE{BOX.END2,ison2} > 0
      fullb(ison2) = 1;
      TREE{BOX.C2,ibox} = TREE{BOX.C2,ibox} - nz; 
   else
       TREE{BOX.C2,ibox} = -1; 
       nz = nz+1; 
   end
end

TREE = TREE(:,fullb==1); 

TREE{BOX.END1,4} = 1; sumb = TREE{BOX.END2,4}; 
for ibox = 5:(7-nz) 
   TREE{BOX.END1,ibox} = sumb+1; 
   sumb = sumb + TREE{BOX.END2,ibox}; 
end

nbox_max = 112;
% parameters to keep subdividing the tree, if necessary
ibox_last = 3;
ncreated  = sum(fullb(4:7));
ilevel    = 2;

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

%nboxes = ibox_last;
nboxes = size(TREE,2);
NODES  = cell(50,nboxes);
NODES(1:7,:) = TREE;

return

