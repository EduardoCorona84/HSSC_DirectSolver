function K_HSS = HSS1D_compress_rectangular_nsym(type,Xin,Xout,TREE,acc,params,dim)
%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%     
%
%     Version of Martinsson's 1D Solver Code. Given points on a few layers on a 
%     curve X, kernel and parameters, it generates an HSS compressed form of
%     the associated matrix K[Xin,Xout]. 
%     
%     K[Xin,Xout] is an M x N rectangular matrix. Different binary trees are used 
%     for Xin and Xout. For simplicity we assume they have the same structure and
%     depth, but this can be generalized. 
%     
%     Inputs: 
%     X  - (N x 2) array with sorted points in a few layers along a curve (usually 
%                  a box boundary or an interface between two boxes). 
%     nbox_max - (int) maximum number of points allowed on a leaf on the tree hierarchy
%     acc      - (double) desired accuracy
%     params   - (struct) parameters
%     dim      - (int) dimension of kernel
%     
%     Outputs: 
%     K_HSS (cell 50 x nboxes) Cell containing HSS tree and skeleton information: 
%         Tree_out: 
%         K_HSS{2,ibox} - level
%         K_HSS{3,ibox} - parent
%         K_HSS{4:5,ibox} - children c1 and c2
%         K_HSS{6:7,ibox} - outgoing index endpoints
%         K_HSS{BOX.NSKEL,ibox}   - nskel 
%         K_HSS{9,ibox}   - kskel
%     
%         Tree_in: 
%         (optional - add level, parent and children info on 12-15)
%         K_HSS{16:17,ibox} - incoming index endpoints
%         K_HSS{18,ibox}   - nskel 
%         K_HSS{19,ibox}   - kskel
%         
%         K_HSS{BOX.I_SKUP,ibox} - I_skup outgoing index with skeleton points
%         K_HSS{BOX.I_SKDN,ibox} - I_skdn incoming index with skeleton points
%     
%         K_HSS{BOX.T_UP,ibox} - Tup outgoing interpolation matrix
%         K_HSS{BOX.J_UP,ibox} - Jup outgoing skeleton permutation vector
%         K_HSS{BOX.T_DN,ibox} - Tdn incoming interpolation matrix
%         K_HSS{BOX.J_DN,ibox} - Jdn incoming skeleton permutation vector
%     
%         K_HSS{BOX.M_SELF,ibox} - Self-Interaction matrices (leaves only) 
%         K_HSS{BOX.M_SIB,ibox} - Sibling Interaction matrices (except top)
%     
%         K_HSS{1,11,21,23:30,24:39,41:45,47:50,ibox} - empty
%


global BOX;       

if nargin < 7
    dim = 1; 
end

% The following parameter can be tuned. 
% Setting it to 1.5 has proved a good balance in many cases.
radius_rel = 1.75;                      
nproxy = 50;       

if strcmp(type,'green')
    mtot   = size(Xin,1);
    ntot   = size(Xout,1); 
else
    K = Xin; 
    [mtot,ntot] = size(K);
    clear Xin Xout; 
end

% Compute the tree structure.
if length(TREE) == 1
    nbox_max = TREE;
    TREE_out  = LOCAL_get_tree(ntot,nbox_max); 
    TREE_in  = LOCAL_get_tree_in(mtot,TREE_out);  
else
    TREE_out = TREE(1:7,:); 
    TREE_in  = TREE(11:17,:);  
end

nboxes = size(TREE_out,2); 
K_HSS = cell(50,nboxes); 
K_HSS(1:7,1:nboxes)  = TREE_out; 
K_HSS(16:17,1:nboxes) = TREE_in(6:7,1:nboxes); 

if strcmp(type,'green')
% Construct the list of neighbors of any cell.
%LIST_NEI_in  = LOCAL_construct_potential_neighborlist(Xin,TREE_in,radius_rel);
%LIST_NEI_out = LOCAL_construct_potential_neighborlist(Xout,TREE_out,radius_rel);
[LIST_NEI_in,LIST_NEI_out] = LOCAL_construct_potential_neighborlist(Xin,Xout,K_HSS,radius_rel);  

% nproxy ~ 4*(2*pi*r)/(2*pi/k) + 50 = 4*k*r + 50, so there are 4 points per wavelength
% (for oscillatory problems).
if strcmp(params.flag_pot,'SL_L_3D')
    nproxy_in = 4*nproxy; 
    nproxy_out = 4*nproxy;       
else
    nproxy_in =  nproxy; 
    nproxy_out = nproxy;    
end 

else
    depth = K_HSS{BOX.LEVEL,nboxes}; 
    I_src = cell(depth+1,2);
end  

% Compress all K_HSS, going from smaller to larger.
for ibox = nboxes:(-1):2

  % Construct the index vectors for ibox.
  if ( (K_HSS{BOX.C1,ibox}<=0) & (K_HSS{BOX.C2,ibox}<=0) ) % ibox has no sons.
    indskel_out = dim*K_HSS{BOX.END1,ibox} - dim + (1:dim*K_HSS{BOX.END2,ibox});
    indskel_in  = dim*K_HSS{16,ibox} - dim + (1:dim*K_HSS{17,ibox});
  else % ibox has two sons
    ison1       = K_HSS{BOX.C1,ibox};
    ison2       = K_HSS{BOX.C2,ibox};
    indskel_out = [K_HSS{BOX.I_SKUP,ison1},K_HSS{BOX.I_SKUP,ison2}];
    indskel_in  = [K_HSS{BOX.I_SKDN,ison1},K_HSS{BOX.I_SKDN,ison2}];
  end
  
  if strcmp(type,'green')
        % Determine a circle that circumscribes the cell:
        [xxc,R] = LOCAL_get_circum_circle([Xin(sort(indskel_in),:) ; Xout(sort(indskel_out),:)]);
        xxc_in = xxc; xxc_out = xxc;   
        R_in = R; R_out = R; 
        
        if strcmp(params.flag_pot,'SL_H_2D') | strcmp(params.flag_pot,'SL_H_3D')
            nproxy_in = min(max(nproxy,round(4*params.kh*radius_rel*R_in + 30)),mtot);
            nproxy_out = min(max(nproxy,round(4*params.kh*radius_rel*R_out + 30)),ntot);
        end
        
        %[xxc_in,R_in] = LOCAL_get_circum_circle(Xin(sort(indskel_in),:));
        %[xxc_out,R_out] = LOCAL_get_circum_circle(Xout(sort(indskel_out),:));
  
        Cproxy_in  = LOCAL_construct_circle(xxc_in, radius_rel*R_in, nproxy_in);   
        Cproxy_out = LOCAL_construct_circle(xxc_out, radius_rel*R_out, nproxy_out);  
  
        % Find all contour points inside the circle:
        indskel_in_maybeinside = [];
        indskel_out_maybeinside = [];
  
        % Build inskel_out_maybeinside
        for jbox = LIST_NEI_in{ibox}
            if ( (K_HSS{BOX.C1,jbox}<=0) & (K_HSS{BOX.C2,jbox}<=0) )
            indskel_out_maybeinside = [indskel_out_maybeinside,dim*(K_HSS{BOX.END1,jbox}) - dim + (1:dim*K_HSS{BOX.END2,jbox})];
            else
            jbox_son1 = K_HSS{BOX.C1,jbox};
            jbox_son2 = K_HSS{BOX.C2,jbox};
            indskel_out_maybeinside = [indskel_out_maybeinside,K_HSS{BOX.I_SKUP,jbox_son1},K_HSS{BOX.I_SKUP,jbox_son2}];
            end
        end
  
        % Build inskel_in_maybeinside
        for jbox = LIST_NEI_out{ibox}  
            if ( (K_HSS{BOX.C1,jbox}<=0) & (K_HSS{BOX.C2,jbox}<=0) )
            indskel_in_maybeinside  = [indskel_in_maybeinside, dim*(K_HSS{16,jbox}) - dim + (1:dim*K_HSS{17,jbox})];
            else
            jbox_son1 = K_HSS{BOX.C1,jbox};
            jbox_son2 = K_HSS{BOX.C2,jbox};
            indskel_in_maybeinside  = [indskel_in_maybeinside, K_HSS{BOX.I_SKDN,jbox_son1},K_HSS{BOX.I_SKDN,jbox_son2}];
            end
        end
  
        relind = find( ((Xout(indskel_out_maybeinside,1)-xxc_in(1)).^2 + ...
                        (Xout(indskel_out_maybeinside,2)-xxc_in(2)).^2) < ((radius_rel*R_in)^2) );
        
        indskel_out_inside = indskel_out_maybeinside(relind);
        
        relind = find( ((Xin(indskel_in_maybeinside,1) -xxc_out(1)).^2 + ...
                        (Xin(indskel_in_maybeinside,2) -xxc_out(2)).^2) < ((radius_rel*R_out)^2) );
        
        indskel_in_inside  = indskel_in_maybeinside(relind);
  
        Xin_proxy = [Cproxy_in(1,:)' Cproxy_in(4,:)'];   
        Xout_proxy = [Cproxy_out(1,:)' Cproxy_out(4,:)'];
  
        params.proxy = 1; 
        A_out = [Kernel_Eval(Xin(indskel_in_inside,:),Xout(indskel_out,:),params);...
              Kernel_Eval(Xin_proxy,Xout(indskel_out,:),params)]; 
        A_in = [Kernel_Eval(Xin(indskel_in,:),Xout(indskel_out_inside,:),params),...
              Kernel_Eval(Xin(indskel_in,:),Xout_proxy,params)]; 
  else
      lev = K_HSS{BOX.LEVEL,ibox}; 
      if ( (K_HSS{BOX.C1,ibox}<=0) && (K_HSS{BOX.C2,ibox}<=0) )
        % Neutered index (I \ I_i)
         ind_offd_out = [1:(K_HSS{BOX.END1,ibox}-1),(K_HSS{BOX.END1,ibox}+K_HSS{BOX.END2,ibox}):ntot];
         ind_offd_in = [1:(K_HSS{BOX.END1+10,ibox}-1),(K_HSS{BOX.END1+10,ibox}+K_HSS{BOX.END2+10,ibox}):mtot]; 
      else
          % Neutered index U_{j \neq i} {I_j^{src}}   
         ind_offd_out = I_src{lev+1,1};
         ind_offd_in  = I_src{lev+1,2}; 
         
         ind_offd_out(ind_offd_out>=K_HSS{BOX.END1,ibox} & ...
                  ind_offd_out<K_HSS{BOX.END1,ibox}+K_HSS{BOX.END2,ibox}) = [];     
         ind_offd_in(ind_offd_in>=K_HSS{BOX.END1+10,ibox} & ...
                  ind_offd_in<K_HSS{BOX.END1+10,ibox}+K_HSS{BOX.END2+10,ibox}) = [];        
      end
        
        % HSS block-column
        A_out = K(ind_offd_in,indskel_out);
        % HSS block-row
        A_in = K(indskel_in,ind_offd_out);
  end
      
  % Compute the skeletons.
  [Tout,Jout] = ID(A_out, acc);
  [Tin, Jin ] = ID(A_in.',acc);
  k_out = size(Tout,1); 
  k_in  = size(Tin,1);
  
  % Record the outgoing skeletons:
  K_HSS{BOX.I_SKUP,ibox} = indskel_out(Jout(1:k_out));
  K_HSS{BOX.T_UP,ibox} = Tout;
  K_HSS{BOX.J_UP,ibox} = Jout;
  % Record the incoming skeletons:
  K_HSS{BOX.I_SKDN,ibox} = indskel_in(Jin(1:k_in));
  K_HSS{BOX.T_DN,ibox} = Tin;
  K_HSS{BOX.J_DN,ibox} = Jin;
  % Record nskel and kskel
  K_HSS{BOX.NSKEL,ibox} = length(indskel_out);
  K_HSS{BOX.KSKEL,ibox} = k_out;
  K_HSS{ 18,ibox} = length(indskel_in);
  K_HSS{ 19,ibox} = k_in;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if strcmp(type,'brute') & K_HSS{BOX.LEVEL,ibox}>0
     lev = K_HSS{BOX.LEVEL,ibox};  
     I_src{lev,1} = [I_src{lev,1} K_HSS{BOX.I_SKUP,ibox} ];     
     I_src{lev,2} = [I_src{lev,2} K_HSS{BOX.I_SKDN,ibox} ];     
  end
end

% Construct the matrices representing self-interactions on the leaves and
% sibling interactions on children boxes
for ibox = nboxes:(-1):1
  if ( (K_HSS{BOX.C1,ibox}<=0) & (K_HSS{BOX.C2,ibox}<=0) )
    ind_out = dim*K_HSS{BOX.END1,ibox} - dim + (1:dim*K_HSS{BOX.END2,ibox});
    ind_in  = dim*K_HSS{16,ibox} - dim + (1:dim*K_HSS{17,ibox});
    params.proxy = 0;
    
    if strcmp(type,'green')
        K_HSS{BOX.M_SELF,ibox} = Kernel_Eval(Xin(ind_in,:),Xout(ind_out,:),params); 
    else
        K_HSS{BOX.M_SELF,ibox} = K(ind_in,ind_out); 
    end
    
  elseif (K_HSS{BOX.C1,ibox}>0) & (K_HSS{BOX.C2,ibox}>0)
    ison1 = K_HSS{BOX.C1,ibox};
    ison2 = K_HSS{BOX.C2,ibox};
    params.proxy = 0; 
    
    if strcmp(type,'green')
        K_HSS{BOX.M_SIB,ison1} = Kernel_Eval(Xin(K_HSS{BOX.I_SKDN,ison1},:),Xout(K_HSS{BOX.I_SKUP,ison2},:),params); 
        K_HSS{BOX.M_SIB,ison2} = Kernel_Eval(Xin(K_HSS{BOX.I_SKDN,ison2},:),Xout(K_HSS{BOX.I_SKUP,ison1},:),params); 
    else
        K_HSS{BOX.M_SIB,ison1} = K(K_HSS{BOX.I_SKDN,ison1},K_HSS{BOX.I_SKUP,ison2}); 
        K_HSS{BOX.M_SIB,ison2} = K(K_HSS{BOX.I_SKDN,ison2},K_HSS{BOX.I_SKUP,ison1});  
    end
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TREE = LOCAL_get_tree(ntot,nbox_max)
global BOX
TREE = cell(7,ntot);

% Construct the top node.
TREE{2,1} = 0;
TREE{3,1} = NaN;
TREE{4,1} = -1;
TREE{5,1} = -1;
TREE{6,1} = 1;
TREE{7,1} = ntot;
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
    nbox = TREE{7,ibox};
    if (nbox > nbox_max)
      nhalf             = ceil(nbox/2);
      ibox_son1         = ibox_last + ncreated + 1;
      ibox_son2         = ibox_last + ncreated + 2;
      TREE{4,ibox}      = ibox_son1;
      TREE{5,ibox}      = ibox_son2;
      TREE{2,ibox_son1} = ilevel+1;
      TREE{2,ibox_son2} = ilevel+1;
      TREE{3,ibox_son1} = ibox;
      TREE{3,ibox_son2} = ibox;
      TREE{4,ibox_son1} = -1;
      TREE{4,ibox_son2} = -1;
      TREE{5,ibox_son1} = -1;
      TREE{5,ibox_son2} = -1;
      TREE{6,ibox_son1} = TREE{6,ibox};
      TREE{6,ibox_son2} = TREE{6,ibox} + nhalf;
      TREE{7,ibox_son1} = nhalf;
      TREE{7,ibox_son2} = TREE{7,ibox} - nhalf;
      ncreated          = ncreated + 2;
    end
  end
  ilevel = ilevel + 1;
end

nboxes = ibox_last;
TREE = TREE(:,1:nboxes); 

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TREE_in = LOCAL_get_tree_in(mtot,TREE_out)
global BOX
nboxes = size(TREE_out,2); 
TREE_in = cell(7,nboxes);
TREE_in(1:5,:) = TREE_out(1:5,:); 

% Construct the top node.
TREE_in{6,1} = 1;
TREE_in{7,1} = mtot;

for ibox = 1:nboxes
   mbox = TREE_in{7,ibox}; 
   if (TREE_out{4,ibox}>0 & TREE_in{5,ibox}>0) % box has two sons
       mhalf = ceil(mbox/2); 
       ison1 = TREE_out{4,ibox};
       ison2 = TREE_out{5,ibox};
       TREE_in{6,ison1} = TREE_in{6,ibox}; 
       TREE_in{6,ison2} = TREE_in{6,ibox} + mhalf; 
       TREE_in{7,ison1} = mhalf; 
       TREE_in{7,ison2} = TREE_in{7,ibox} - mhalf; 
   elseif (TREE_out{4,ibox}>0 | TREE_in{5,ibox}>0) % box has one son
       ison = max(TREE_out{4,ibox},TREE_out{5,ibox});
       TREE_in{6,ison} = TREE_in{6,ibox}; 
       TREE_in{7,ison} = TREE_in{7,ison}; 
   end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xxc, R] = LOCAL_get_circum_circle(X)

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

function C = LOCAL_construct_circle(xxc, R, n)

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

function [LIST_NEI_in,LIST_NEI_out] = LOCAL_construct_potential_neighborlist(X_in,X_out,K_HSS,radius_rel)
global BOX
nboxes = size(K_HSS,2);

% Note that in the construction, we temporarily include
% a box in its own list of neighbors.
LIST_NEI_in    = cell(1,nboxes);
LIST_NEI_in{1} = [1];
LIST_NEI_out = LIST_NEI_in; 

for ibox = 2:nboxes
    ind_in  = K_HSS{16,ibox}-1+(1:K_HSS{17,ibox}); 
    ind_out =  K_HSS{BOX.END1,ibox}-1+(1:K_HSS{BOX.END2,ibox}); 
  %[xxc_in, R_in] = LOCAL_get_circum_circle(X_in(ind_in,:));
  %[xxc_out, R_out] = LOCAL_get_circum_circle(X_out(ind_out,:)); 
  [xxc,R] = LOCAL_get_circum_circle([X_in(ind_in,:) ; X_out(ind_out,:)]);
  xxc_in = xxc; xxc_out = xxc; R_in = R; R_out = R;          
  
  
  % Collect of list of all the father's neighbors' sons:
  potneis_in = []; potneis_out = [];
  
  for jbox = LIST_NEI_in{K_HSS{3,ibox}};
    if (K_HSS{BOX.C1,jbox} > 0)
      potneis_in = [potneis_in,K_HSS{BOX.C1,jbox}];
    end
    if (K_HSS{BOX.C2,jbox} > 0)
      potneis_in = [potneis_in,K_HSS{BOX.C2,jbox}];
    end
  end
  
  for jbox = LIST_NEI_out{K_HSS{3,ibox}};
    if (K_HSS{BOX.C1,jbox} > 0)
      potneis_out = [potneis_out,K_HSS{BOX.C1,jbox}];
    end
    if (K_HSS{BOX.C2,jbox} > 0)
      potneis_out = [potneis_out,K_HSS{BOX.C2,jbox}];
    end
  end
  
  % Loop over all potential neighbors to check which ones are actually "close".
  for jbox = potneis_in
    ind    = K_HSS{BOX.END1,jbox} - 1 + (1:K_HSS{BOX.END2,jbox});
    distsq = (X_out(ind,1)-xxc_in(1)).^2 + (X_out(ind,2)-xxc_in(2)).^2;
    if (min(distsq) < radius_rel*radius_rel*R_in*R_in)
      LIST_NEI_in{ibox} = [LIST_NEI_in{ibox},jbox];
    end
  end
  
  for jbox = potneis_out
    ind    = K_HSS{16,jbox} - 1 + (1:K_HSS{17,jbox});
    distsq = (X_in(ind,1)-xxc_out(1)).^2 + (X_in(ind,2)-xxc_out(2)).^2;
    if (min(distsq) < radius_rel*radius_rel*R_out*R_out)   
      LIST_NEI_out{ibox} = [LIST_NEI_out{ibox},jbox];
    end
  end
end
% Finally we remove the box itself from the list of neighbors.
for ibox = 2:nboxes
  j = find(LIST_NEI_in{ibox} == ibox);
  LIST_NEI_in{ibox} = sort(LIST_NEI_in{ibox}([1:(j-1),(j+1):end]));
  
  j = find(LIST_NEI_out{ibox} == ibox);
  LIST_NEI_out{ibox} = sort(LIST_NEI_out{ibox}([1:(j-1),(j+1):end]));   
end

return