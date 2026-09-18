function K_HSS = HSS1D_compress_submatrix_fsym_giventree(X,TREE,Box_dn,Box_up,acc,params,dim)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        K_HSS = HSS1D_compress_submatrix_fsym_giventree(X,TREE,Box_dn,Box_up,acc,params,dim)
%
%    DESCRIPTION:
%        Given points on a few layers on a curve X, symmetric kernel and parameters, it generates
%        an HSS compressed form of a submatrix of K[X,X] given the tree structure in TREE.
%
%    INFILE VARIABLES:
%        radius_rel (double) specifies the relative radius of the proxy circles. default = 1.75.
%        nproxy (int) specifies the number of samples to be used on the proxy circles. default = 50.
%
%    INPUT:
%        X is an (N x 2) array with sorted points in a few layers along a curve (usually a box
%            boundary or an interface between two boxes).
%        TREE is a {50 x nboxes} cell array that gives the tree structure. The submatrix that is
%            compressed corresponds to TREE(:,Box_dn) and TREE(:,Box_up).
%        Box_dn / Box_up are the (nboxes x 1) boolean arrays that determine the row / column
%            subtrees.
%        acc (double) is the desired accuracy.
%        params is a struct of kernel evaluation parameters that can be generated using
%            ../Common/HSS_tree_parameters.m.
%        dim (int) is the dimension of kernel (scalar -> dim = 1)
%
%    OUTPUT:
%        K_HSS is the cell array containing the HSS tree and skeleton information
%            (see ../Common/box_constants for details)
%


global BOX 

if nargin < 7
    dim = 1; 
end

% The following parameter can be tuned. 
% Setting it to 1.5 has proved a good balance in many cases.
radius_rel = 1.75;
nproxy = 50;
ntot   = size(X,1);

% Compute the tree strcture.
K_HSS  = TREE; 
nboxes = size(K_HSS,2);

% Construct the list of neighbors of any cell.
LIST_NEI = LOCAL_construct_potential_neighborlist(X,K_HSS,radius_rel);

% Compress all K_HSS, going from smaller to larger.
for ibox = nboxes:(-1):2
    if (Box_dn(ibox)==true || Box_up(ibox)==true) 
        % Construct the index vectors for ibox.
        if ( (K_HSS{BOX.C1,ibox}<=0) && (K_HSS{BOX.C2,ibox}<=0) ) % ibox has no sons.
            indskel = dim*K_HSS{BOX.END1,ibox} - dim + (1:dim*K_HSS{BOX.END2,ibox});
            % else, ibox has at least one son
        elseif ( (K_HSS{BOX.C1,ibox}>0) && (K_HSS{BOX.C2,ibox}>0) )
            ison1       = K_HSS{BOX.C1,ibox};
            ison2       = K_HSS{BOX.C2,ibox};
            indskel = [K_HSS{BOX.I_SKUP,ison1},K_HSS{BOX.I_SKUP,ison2}];
        else
            ison = max([K_HSS{BOX.C1,ibox} K_HSS{BOX.C2,ibox}]);
            indskel = K_HSS{BOX.I_SKUP,ison};
        end
  
        % Determine a circle that circumscribes the cell:
        [xxc,R] = LOCAL_get_circum_circle(X(sort(indskel),:));
        
        % nproxy ~ 4*(2*pi*r)/(2*pi/k) + 50 = 4*k*r + 50, so there are 4 points per wavelength
        % (for oscillatory problems).
        if strcmp(params.flag_pot,'SL_H_2D') || strcmp(params.flag_pot,'SL_H_3D')
            nproxy = min(max(50,round(4*params.kh*radius_rel*R + 30)),ntot);  
        end 
        
        Cproxy  = LOCAL_construct_circle(xxc, radius_rel*R, nproxy);
        % Find all contour points inside the circle:
        indskel_maybeinside = [];
        for jbox = LIST_NEI{ibox}
            if ( (K_HSS{BOX.C1,jbox}<=0) && (K_HSS{BOX.C2,jbox}<=0) )
            indskel_maybeinside = [indskel_maybeinside,dim*(K_HSS{BOX.END1,jbox}) - dim + (1:dim*K_HSS{BOX.END2,jbox})];
            else
            jbox_son1 = K_HSS{BOX.C1,jbox};
            jbox_son2 = K_HSS{BOX.C2,jbox};
            
            if jbox_son1>0 && jbox_son2>0
                indskel_maybeinside = [indskel_maybeinside,K_HSS{BOX.I_SKUP,jbox_son1},K_HSS{BOX.I_SKUP,jbox_son2}];
            else
                jbox_son = max([jbox_son1 jbox_son2]); 
                indskel_maybeinside = [indskel_maybeinside,K_HSS{BOX.I_SKUP,jbox_son}];
            end
            end
        end
  
        relind = find( ((X(indskel_maybeinside,1)-xxc(1)).^2 + ...
                      (X(indskel_maybeinside,2)-xxc(2)).^2) < ((radius_rel*R)^2) );
        indskel_inside = indskel_maybeinside(relind);
  
        X_proxy = [Cproxy(1,:)' Cproxy(4,:)']; 
        
        if params.dim == 3
            X_proxy(:,3) = fnval(params.pp,X_proxy.').'; 
        end
  
        A21proxy = [Kernel_Eval(X(indskel_inside,:),X(indskel,:),params);...
                  Kernel_Eval(X_proxy,X(indskel,:),params)]; 
      
        % Compute the skeletons.
        [T,J] = ID(A21proxy, acc);
        k = size(T,1);
  
        % Record the outgoing skeletons:
        K_HSS{BOX.I_SKUP,ibox} = indskel(J(1:k));
        K_HSS{BOX.T_UP,ibox} = T;
        K_HSS{BOX.J_UP,ibox} = J;

        % Record nskel and kskel
        K_HSS{BOX.NSKEL,ibox} = length(indskel);
        K_HSS{BOX.KSKEL,ibox} = k;
    end
end

% Construct the matrices representing self-interactions on the leaves.
for ibox = nboxes:(-1):2
    if ( (K_HSS{BOX.C1,ibox}<=0) && (K_HSS{BOX.C2,ibox}<=0) )
        ind = dim*K_HSS{BOX.END1,ibox} - dim + (1:dim*K_HSS{BOX.END2,ibox});
    
        if (Box_dn(ibox)==true && Box_up(ibox)==true)
            K_HSS{BOX.M_SELF,ibox} = Kernel_Eval(X(ind,:),X(ind,:),params);
        else
            K_HSS{BOX.M_SELF,ibox} = sparse(K_HSS{BOX.END2,ibox},K_HSS{BOX.END2,ibox});
        end
    end
end

% Construct the matrices for sibling interactions.
for ibox = nboxes:(-1):1
  if ( (K_HSS{BOX.C1,ibox}>0) && (K_HSS{BOX.C2,ibox}>0) )
    ison1 = K_HSS{BOX.C1,ibox};
    ison2 = K_HSS{BOX.C2,ibox};
    
    if ((Box_dn(ison1)==true && Box_up(ison2)==true) || (Box_dn(ison2)==true && Box_up(ison1)==true))
        
        K_HSS{BOX.M_SIB,ison1} = Kernel_Eval(X(K_HSS{BOX.I_SKUP,ison1},:),X(K_HSS{BOX.I_SKUP,ison2},:),params);
    else
        K_HSS{BOX.M_SIB,ison1} = sparse(K_HSS{BOX.KSKEL,ison1},K_HSS{BOX.KSKEL,ison2});
    end
    
    K_HSS{BOX.M_SIB,ison2} = K_HSS{BOX.M_SIB,ison1}.';
    
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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