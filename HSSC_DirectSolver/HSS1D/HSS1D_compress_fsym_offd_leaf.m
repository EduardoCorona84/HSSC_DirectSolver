function KHSS = HSS1D_compress_fsym_offd_leaf(X,TREE,B1,acc,params,dim)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        KHSS = HSS1D_compress_fsym_offd_leaf(X,TREE,B1,acc,params,dim)
%
%    DESCRIPTION:
%        This function compresses the symmetric kernel function matrix K[X,X] into HSS form (K_HSS)
%        respective of the given parameters. This function produces a matrix that keeps only off-
%        diagonal blocks for each box at leaf level, determined by parameter B1.
%
%    INFILE VARIABLES:
%        radius_rel (double) specifies the relative radius of the proxy circles. default = 1.75.
%        nproxy (int) specifies the number of samples to be used on the proxy circles. default = 50.
%
%    INPUT:
%        X is an (N x 2) array with sorted points in a few layers along a curve (usually a box
%            boundary or an interface between two boxes).
%        TREE is a {50 x nboxes} cell array that gives the underlying binary tree structure of KHSS.
%        B1 is a (N x 1) boolean array where for each box i, letting B1:=B1(I_src(i)), the self
%            interaction matrices are
%                          D = [0  D(B1i,~B1i]
%                              [D(~B1i,B1i) 0]
%            Interpolation operators R and L are also zeroed out accordingly.
%        acc (double) is the desired accuracy.
%        params is a struct of kernel evaluation parameters that can be generated using
%            ../Common/HSS_tree_parameters.m
%        dim (int) is the dimension of the kernel (scalar -> dim = 1)
%
%    OUTPUT:
%        KHSS is the cell array containing the HSS tree and skeleton information
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
ntot   = size(X,1);

% Compute the tree strcture.
KHSS  = TREE; 
nboxes = size(KHSS,2);

% Construct the list of neighbors of any cell.
LIST_NEI = LOCAL_construct_potential_neighborlist(X,KHSS,radius_rel);

% Compress all KHSS, going from smaller to larger.
for ibox = nboxes:(-1):2

  % Construct the index vectors for ibox.
  if ( (KHSS{BOX.C1,ibox}<=0) && (KHSS{BOX.C2,ibox}<=0) ) % ibox has no sons.
      indskel = dim*KHSS{BOX.END1,ibox} - dim + (1:dim*KHSS{BOX.END2,ibox});
      % B1_src
      KHSS{23,ibox}  = B1(indskel); 
    % else, ibox has at least one son
  elseif ( (KHSS{BOX.C1,ibox}>0) && (KHSS{BOX.C2,ibox}>0) )
      ison1       = KHSS{BOX.C1,ibox};
      ison2       = KHSS{BOX.C2,ibox};
      indskel = [KHSS{BOX.I_SKUP,ison1},KHSS{BOX.I_SKUP,ison2}];
      % B1_src
      %KHSS{23,ibox}  = [KHSS{24,ison1},KHSS{24,ison2}];
  else
      ison = max([KHSS{BOX.C1,ibox} KHSS{BOX.C2,ibox}]);
      indskel = KHSS{BOX.I_SKUP,ison};
      % B1_src
      %KHSS{23,ibox} = KHSS{24,ison}; 
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
    if ( (KHSS{BOX.C1,jbox}<=0) && (KHSS{BOX.C2,jbox}<=0) )
      indskel_maybeinside = [indskel_maybeinside,dim*(KHSS{BOX.END1,jbox}) - dim + (1:dim*KHSS{BOX.END2,jbox})];
    else
      jbox_son1 = KHSS{BOX.C1,jbox};
      jbox_son2 = KHSS{BOX.C2,jbox};
      if (jbox_son1>0 && jbox_son2>0)
          indskel_maybeinside = [indskel_maybeinside,KHSS{BOX.I_SKUP,jbox_son1},KHSS{BOX.I_SKUP,jbox_son2}];
      else
          jbox_son = max(jbox_son1,jbox_son2);
          indskel_maybeinside = [indskel_maybeinside,KHSS{BOX.I_SKUP,jbox_son}];
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
  KHSS{BOX.I_SKUP,ibox} = indskel(J(1:k));
  %KHSS{24,ibox} = KHSS{23,ibox}(J(1:k)); 
  KHSS{BOX.T_UP,ibox} = T;
  KHSS{BOX.J_UP,ibox} = J;

  % Record nskel and kskel
  KHSS{BOX.NSKEL,ibox} = length(indskel);
  KHSS{BOX.KSKEL,ibox} = k;
end

%Self and Sibling Interactions
for ibox = nboxes:(-1):1
  if ( (KHSS{BOX.C1,ibox}<=0) && (KHSS{BOX.C2,ibox}<=0) )
    ind = dim*KHSS{BOX.END1,ibox} - dim + (1:dim*KHSS{BOX.END2,ibox});
    B1_src = KHSS{23,ibox}; B2_src = ~B1_src; 
    
    KHSS{BOX.M_SELF,ibox} = zeros(KHSS{BOX.END2,ibox},KHSS{BOX.END2,ibox}); 
    if min(sum(B1_src),sum(B2_src))>0
    KHSS{BOX.M_SELF,ibox}(B1_src,B2_src) = Kernel_Eval(X(ind(B1_src),:),X(ind(B2_src),:),params); 
    KHSS{BOX.M_SELF,ibox}(B2_src,B1_src) = KHSS{BOX.M_SELF,ibox}(B1_src,B2_src).'; 
    end
  else
    ison1 = KHSS{BOX.C1,ibox};
    ison2 = KHSS{BOX.C2,ibox};
    if ison1>0 && ison2>0 
       KHSS{BOX.M_SIB,ison1} = Kernel_Eval(X(KHSS{BOX.I_SKUP,ison1},:),X(KHSS{BOX.I_SKUP,ison2},:),params);
       KHSS{BOX.M_SIB,ison2} = KHSS{BOX.M_SIB,ison1}.'; 
    end
  end
end

% Have to modify sibling interactions and interpolation operators to
% isolate off-diagonal part. This step should be eventually done directly. 
for ibox = nboxes:(-1):1
    ison1 = KHSS{BOX.C1,ibox};
    ison2 = KHSS{BOX.C2,ibox};
    
    if ibox>1
    R = [eye(KHSS{BOX.KSKEL,ibox}) KHSS{BOX.T_UP,ibox}]; R(:,KHSS{BOX.J_UP,ibox}) = R;  
    
    if ( ison1<=0 && ison2<=0 )
      indskel = dim*KHSS{BOX.END1,ibox} - dim + (1:dim*KHSS{BOX.END2,ibox});
      n1 = sum(KHSS{23,ibox}); n2 = KHSS{BOX.NSKEL,ibox} - n1;
      Rhat = [R(:,1:n1) zeros(KHSS{BOX.KSKEL,ibox},n2) ; zeros(KHSS{BOX.KSKEL,ibox},n1) R(:,n1+1:end)]; 
    else
      if ison1>0 && ison2>0
          indskel = [KHSS{BOX.I_SKUP,ison1},KHSS{BOX.I_SKUP,ison2}];
          R_c1 = KHSS{BOX.RHO,ison1}; R_c2 = KHSS{BOX.RHO,ison2}; 
          k1 = size(R_c1,1)/2; k2 = size(R_c2,1)/2; 
          
          Rhat = [R(:,1:k1)*R_c1(1:k1,:) R(:,k1+1:end)*R_c2(1:k2,:) ; ...
                  R(:,1:k1)*R_c1(k1+1:end,:) R(:,k1+1:end)*R_c2(k2+1:end,:)];
      else
          ison = max(ison1,ison2); 
          indskel = KHSS{BOX.I_SKUP,ison}; 
          R_c = KHSS{BOX.RHO,ison}; k = size(R_c,1)/2; 
          
          Rhat = [R*R_c(1:k,:) ; R*R_c(k+1:end,:) ]; 
      end
    end
    
    [T,J] = ID(Rhat,acc); k = size(T,1); 
    KHSS{BOX.I_SKUP,ibox} = indskel(J(1:k)); 
    KHSS{BOX.T_UP,ibox} = T; 
    KHSS{BOX.J_UP,ibox} = J; 
    KHSS{BOX.RHO,ibox} = Rhat(:,J(1:k));
    KHSS{BOX.NSKEL,ibox} = length(indskel); 
    KHSS{BOX.KSKEL,ibox} = k; 
    end
    
    if ison1>0 && ison2>0 
        k1 = size(KHSS{BOX.M_SIB,ison1},1); k2 = size(KHSS{BOX.M_SIB,ison1},2); 
        KHSS{BOX.M_SIB,ison1} = KHSS{BOX.RHO,ison1}(1:k1,:).'*KHSS{BOX.M_SIB,ison1}*KHSS{BOX.RHO,ison2}(k2+1:end,:) + ...
                        KHSS{BOX.RHO,ison1}(k1+1:end,:).'*KHSS{BOX.M_SIB,ison1}*KHSS{BOX.RHO,ison2}(1:k2,:);
        KHSS{BOX.M_SIB,ison2} = KHSS{BOX.M_SIB,ison1}.';             
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

function LIST_NEI = LOCAL_construct_potential_neighborlist(X,KHSS,radius_rel)
global BOX 

nboxes = size(KHSS,2);

% Note that in the construction, we temporarily include
% a box in its own list of neighbors.
LIST_NEI    = cell(1,nboxes);
LIST_NEI{1} = [1];
for ibox = 2:nboxes
  [xxc, R] = LOCAL_get_circum_circle(X(KHSS{BOX.END1,ibox}-1+(1:KHSS{BOX.END2,ibox}),:));
  % Collect of list of all the father's neighbors' sons:
  potneis = [];
  for jbox = LIST_NEI{KHSS{3,ibox}};
    if (KHSS{BOX.C1,jbox} > 0)
      potneis = [potneis,KHSS{BOX.C1,jbox}];
    end
    if (KHSS{BOX.C2,jbox} > 0)
      potneis = [potneis,KHSS{BOX.C2,jbox}];
    end
  end
  % Loop over all potential neighbors to check which ones are actually "close".
  for jbox = potneis
    ind    = KHSS{BOX.END1,jbox} - 1 + (1:KHSS{BOX.END2,jbox});
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