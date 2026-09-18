function K_HSS = HSS1D_compress_green_nsym_Kproxy(Xin,Xout,TREE,acc,params,dim)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        K_HSS = HSS1D_compress_green_nsym_Kproxy(Xin,Xout,TREE,acc,params,dim)
%
%    DESCRIPTION:
%        This function generates an HSS compressed form of the associated kernel matrix K[Xin,Xout].
%        K[Xin,Xout] is assumed to be a square matrix, where Xin is not necessarily equal to Xout.
%        A uniform binary tree is used for both sets of points.
%
%    INFILE VARIABLES:
%        radius_rel is an infile variable that can be tuned. It specifies the relative radius of
%            the proxy circles. default = 1.75;
%        nproxy is an infile variable that can be tuned. It specifies the number of samples to be
%            used on the proxy circles. default = 50;
%
%    INPUT:
%        Xin is a (N x 2) matrix of the positions of the points along the curve "acting" on Xout
%        Xout is a (N x 2) matrix of the positions of the points along the curve being "acted upon."
%        TREE can be a {7 x nboxes} cell array that gives the underlying binary tree structure of
%            K_HSS.
%        acc (double) is the desired accuracy.
%        params is a struct of kernel evaluation parameters that can be generated using
%            ../Common/HSS_tree_parameters.m
%        dim (int) is the dimension of the kernel (scalar -> dim = 1)
%
%    OUTPUT:
%        K_HSS is the cell array containing the HSS tree and skeleton information
%            (see ../Common/box_constants for details)
%


global BOX

if nargin < 6
    dim = 1; 
end

% The following parameter can be tuned. 
% Setting it to 1.5 has proved a good balance in many cases.
radius_rel = 1.75;
nproxy = 50;

% Compute the tree strcture.
ntot   = size(Xin,1);

if length(TREE) == 1
    nbox_max = TREE; 
    K_HSS  = LOCAL_get_tree(ntot,nbox_max);
    nboxes = size(K_HSS,2);
else 
    nboxes = size(TREE,2); 
    K_HSS  = cell(50,nboxes);
    K_HSS(BOX.DATA,:) = TREE(BOX.DATA,:); 
    clear TREE;
end

% Construct the list of neighbors of any cell.
LIST_NEI = LOCAL_construct_potential_neighborlist(Xin,K_HSS,radius_rel);

% Compress all K_HSS, going from smaller to larger.
for ibox = nboxes:(-1):2

  % Construct the index vectors for ibox.
  if ( (K_HSS{BOX.C1,ibox}<=0) && (K_HSS{BOX.C2,ibox}<=0) ) % ibox has no sons.
      ind            = dim*K_HSS{BOX.END1,ibox} - dim + (1:dim*K_HSS{BOX.END2,ibox});
    indskel_out = ind;
    indskel_in  = ind;
  else % ibox has two sons
    ison1       = K_HSS{BOX.C1,ibox};
    ison2       = K_HSS{BOX.C2,ibox};
    indskel_out = [K_HSS{BOX.I_SKUP,ison1},K_HSS{BOX.I_SKUP,ison2}];
    indskel_in  = [K_HSS{BOX.I_SKDN,ison1},K_HSS{BOX.I_SKDN,ison2}];
  end
  
  % Determine a circle that circumscribes the cell:
  [xxc,R] = LOCAL_get_circum_circle([Xin(sort(indskel_in),:) ; Xout(sort(indskel_out),:)]);
  
  % nproxy ~ 4*(2*pi*r)/(2*pi/k) + 50 = 4*k*r + 50, so there are 4 points per wavelength
  % (for oscillatory problems).
  if strcmp(params.flag_pot,'SL_H_2D') | strcmp(params.flag_pot,'SL_H_3D')
      nproxy_in = min(max(nproxy,round(4*params.kh*radius_rel*R + 30)),mtot);
      nproxy_out = min(max(nproxy,round(4*params.kh*radius_rel*R + 30)),ntot);
  elseif strcmp(params.flag_pot,'SL_L_3D')
      nproxy_in = 4*nproxy; 
      nproxy_out = 4*nproxy;       
  else
      nproxy_in =  nproxy; 
      nproxy_out = nproxy;    
  end    
  
  Cproxy_in  = LOCAL_construct_circle(xxc, radius_rel*R, nproxy_in);
  Cproxy_out = LOCAL_construct_circle(xxc, radius_rel*R, nproxy_out);
  
  % Find all contour points inside the circle:
  indskel_in_maybeinside = [];
  indskel_out_maybeinside = [];
  for jbox = LIST_NEI{ibox}
    if ( (K_HSS{BOX.C1,jbox}<=0) && (K_HSS{BOX.C2,jbox}<=0) )
      indskel_out_maybeinside = [indskel_out_maybeinside,dim*(K_HSS{BOX.END1,jbox}) - dim + (1:dim*K_HSS{BOX.END2,jbox})];
      indskel_in_maybeinside  = [indskel_in_maybeinside, dim*(K_HSS{BOX.END1,jbox}) - dim + (1:dim*K_HSS{BOX.END2,jbox})];
    else
      jbox_son1 = K_HSS{BOX.C1,jbox};
      jbox_son2 = K_HSS{BOX.C2,jbox};
      indskel_out_maybeinside = [indskel_out_maybeinside,K_HSS{BOX.I_SKUP,jbox_son1},K_HSS{BOX.I_SKUP,jbox_son2}];
      indskel_in_maybeinside  = [indskel_in_maybeinside, K_HSS{BOX.I_SKDN,jbox_son1},K_HSS{BOX.I_SKDN,jbox_son2}];
    end
  end
  
  relind = find( ((Xout(indskel_out_maybeinside,1)-xxc(1)).^2 + ...
                  (Xout(indskel_out_maybeinside,2)-xxc(2)).^2) < ((radius_rel*R)^2) );
  indskel_out_inside = indskel_out_maybeinside(relind);
  relind = find( ((Xin(indskel_in_maybeinside,1) -xxc(1)).^2 + ...
                  (Xin(indskel_in_maybeinside,2) -xxc(2)).^2) < ((radius_rel*R)^2) );
  indskel_in_inside  = indskel_in_maybeinside(relind);

  indtmp = dim*K_HSS{BOX.END1,ibox} - dim + (1:dim*K_HSS{BOX.END2,ibox});
  
  Xin_proxy = [Cproxy_in(1,:).' Cproxy_in(4,:).']; 
  Xout_proxy = [Cproxy_out(1,:)' Cproxy_out(4,:)'];
  
  A21proxy = [Kernel_Eval(Xin(indskel_in_inside,:),Xout(indskel_out,:),params);...
              Kernel_Eval(Xin_proxy,Xout(indskel_out,:),params)]; 
  A12proxy = [Kernel_Eval(Xin(indskel_in,:),Xout(indskel_out_inside,:),params),...
              Kernel_Eval(Xin(indskel_in,:),Xout_proxy,params)]; 
      
  % Compute the skeletons.
  [Tout,Jout] = ID(A21proxy, acc);
  [Tin, Jin ] = ID(A12proxy.',acc);
  k = max(size(Tout,1),size(Tin,1));
  % We sometimes need to enforce that the outgoing rank = incoming rank.
  % (Note that this is implemented rather clumsily.)
  if ~(k == size(Tout,1))
    [Tout,Jout] = ID(A21proxy, k);
  elseif ~(k == size(Tin,1))
    [Tin, Jin ] = ID(A12proxy.',k);
  end
  % Record the outgoing skeletons:
  K_HSS{BOX.I_SKUP,ibox} = indskel_out(Jout(1:k));
  K_HSS{BOX.T_UP,ibox} = Tout;
  K_HSS{BOX.J_UP,ibox} = Jout;
  % Record the incoming skeletons:
  K_HSS{BOX.I_SKDN,ibox} = indskel_in(Jin(1:k));
  K_HSS{BOX.T_DN,ibox} = Tin;
  K_HSS{BOX.J_DN,ibox} = Jin;
  % Record nskel and kskel
  K_HSS{BOX.NSKEL,ibox} = length(indskel_out);
  K_HSS{BOX.KSKEL,ibox} = k;
end

% Construct the matrices representing self-interactions on the leaves.
for ibox = nboxes:(-1):2
  if ( (K_HSS{BOX.C1,ibox}<=0) && (K_HSS{BOX.C2,ibox}<=0) )
    ind = dim*K_HSS{BOX.END1,ibox} - dim + (1:dim*K_HSS{BOX.END2,ibox});
    
    K_HSS{BOX.M_SELF,ibox} = Kernel_Eval(Xin(ind,:),Xout(ind,:),params); 
    %OMNICONT_construct_A_diag(C,ind,flag_pot,params);
  end
end

% Construct the matrices for sibling interactions.
for ibox = nboxes:(-1):1
  if ( (K_HSS{BOX.C1,ibox}>0) && (K_HSS{BOX.C2,ibox}>0) )
    ison1 = K_HSS{BOX.C1,ibox};
    ison2 = K_HSS{BOX.C2,ibox};
    
    K_HSS{BOX.M_SIB,ison1} = Kernel_Eval(Xin(K_HSS{BOX.I_SKDN,ison1},:),Xout(K_HSS{BOX.I_SKUP,ison2},:),params); 
    K_HSS{BOX.M_SIB,ison2} = Kernel_Eval(Xin(K_HSS{BOX.I_SKDN,ison2},:),Xout(K_HSS{BOX.I_SKUP,ison1},:),params); 
    
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xxc, R] = LOCAL_get_circum_circle(X);

nloc       = size(X,1);

% First we TRY to create the circle based on the endpoints.
xxc        = 0.5*(X(1,:) + X(end,:));
distsquare = (X(:,1)-xxc(1)).^2 + (X(:,2)-xxc(2)).^2;
R          = sqrt(max(distsquare));

% The result is an absurdly large circle, then we instead
% base it on the center of mass of the points.
if ( (1.2*R*R) > ( (X(1,1) - X(end,1))^2 + (X(1,2) - X(end,2))^2) )
  xxc        = (1/nloc)*[sum(X(:,1)); sum(X(:,2))];
  distsquare = (X(:,1)-xxc(1)).^2 + (X(:,2)-xxc(2)).^2;
  R          = sqrt(max(distsquare));
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C = LOCAL_construct_circle(xxc, R, n);

tt = linspace(0, 2*pi*(1-1/n), n);
C  = [xxc(1) + R*cos(tt);...
             - R*sin(tt);...
             - R*cos(tt);...
      xxc(2) + R*sin(tt);...
             + R*cos(tt);...
             - R*sin(tt)];

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LIST_NEI = LOCAL_construct_potential_neighborlist(X,K_HSS,radius_rel)
global BOX 

nboxes = size(K_HSS,2);

% Note that in the construction, we temporarily include
% a box in its own list of neighbors.
LIST_NEI    = cell(1,nboxes);
LIST_NEI{1} = [1];
for ibox = 2:nboxes
  [xxc, R] = LOCAL_get_circum_circle(X(K_HSS{BOX.END1,ibox}-1+(1:K_HSS{BOX.END2,ibox}),:));
  % Collect of list of all the father's neighbors' sons:
  potneis = [];
  for jbox = LIST_NEI{K_HSS{3,ibox}};
    if (K_HSS{BOX.C1,jbox} > 0)
      potneis = [potneis,K_HSS{BOX.C1,jbox}];
    end
    if (K_HSS{BOX.C2,jbox} > 0)
      potneis = [potneis,K_HSS{BOX.C2,jbox}];
    end
  end
  % Loop over all potential neighbors to check which ones are actually "close".
  for jbox = potneis
    ind    = K_HSS{BOX.END1,jbox} - 1 + (1:K_HSS{BOX.END2,jbox});
    distsq = (X(ind,1)-xxc(1)).^2 + (X(ind,2)-xxc(2)).^2;
    if (min(distsq) < radius_rel*radius_rel*R*R)
      LIST_NEI{ibox} = [LIST_NEI{ibox},jbox];
    end
  end
end
% Finally we remove the box itself from the list of neighbors.
for ibox = 2:nboxes
  j = find(LIST_NEI{ibox} == ibox);
  LIST_NEI{ibox} = sort(LIST_NEI{ibox}([1:(j-1),(j+1):end]));
end

return