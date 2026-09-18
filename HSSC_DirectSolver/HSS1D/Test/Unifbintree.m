function TREE = Unifbintree(params)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        TREE = Unifbintree(params)
%
%    DESCRIPTION:
%        This file is located in the Test/ folder. The function returns the necessary arrays for a
%        uniformly refined, binary tree (essentially a quadtree) of points in the box [-1,1]^2, of
%        a fixed depth. It alternates between horizontal and vertical bisection of resulting boxes.
%        Parent and children lists, as well other index sets needed for the HSS solver are stored.
%
%    INPUT:
%        params.depth - int - maximum depth of the quadtree
%        params.X_source - (ns x 2) float array - source points
%        params.h - double - grid spacing
%        params.kh - double - kernel parameter (wave number)
%        params.max_particles - int - maximum number of points in a leaf
%        params.layers - int - number of layers needed for skeletons
%
%    OUTPUT:
%        TREE - struct - containing the binary tree:
%            .numlev - (depth x 1 int array - boxes per level
%            .leaves - (numleaf x 2) int array -  leaves box numbers and level
%            .BOX (per box info)
%                .cent  - (total_box x 2) float array - box center
%                .levbox - (total_box x 1) int array - box level
%                .parent - (total_box x 1) int array 0 box parents (0 for root node)
%                .child - (total_par x 2) int array - box's children (non leaf boxes)
%                .I_src - (nskel x 1) int array - index numbers for points on box before merge
%                .I_sk - (k_skel x 1) int array - index numbers for points on box boundary layers
%                    after merging children (so, minus interface)
%                .n_skel - int - size of I_src
%                .k_skel - int - size of I_sk
%                .J - (nskel x 1) int array - [I_sk, I_rs] index with upwards skeleton first,
%                    residual points second. I_sk and I_rs are sorted so that points are on an
%                    oriented curve of layers.
%


global BOX 

depth = params.depth; 
X_source = params.X_source; 
N = size(X_source,1);

% Array Initialization
nboxes = 2^(depth+1)-1;
cent = zeros(nboxes,2); 
cent_int = cent; 
C = cent; P = cent(:,1); 
levbox = P;

%--------------------------------------------------------------------------
% ROOT NODE

% No parents for root
    P(1) = 0;
% 2 potential children, none so far
    nchild=2;
    C(1,1:nchild) = zeros(1,nchild);
% Center of root at [0,0] (of the square [-1,1]^2)
    cent(1,:) = [0 0];
    cent_int(1,:) = [1 1]; 
    
 % Number of boxes per level
    numlev = zeros(1,depth+1);
    numlev(1)=1;
% Total Number of nodes on the tree
    numnodes = 1;
% Box levels 
    levbox(1) = 0;
% Array that stores box numbers according to level in order to traverse the
% tree
    box_numbers = sparse(depth+1,2^depth);
    box_numbers(1,1) = 1;

% Source Points


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% (I) BINARY TREE CONSTRUCTION

% We traverse the tree while the level of the boxes considered is not
% greater than the maximum assigned depth, and as long as bisection is
% required. 

% On the first pass, the necessary arrays are created for each box going 
% down the binary tree. 

lev=0;
num_leaf=0;
num_bisec = zeros(1,depth+1);
Sindex = cell(1,depth+1);

for i=1:depth+1
    if mod(i-1,2) == 0
        Sindex{i} = sparse(2^((i-1)/2),2^((i-1)/2));
    else 
        Sindex{i} = sparse(2^(i/2),2^(i/2-1));
    end
end

% Downwards pass bisection / tree build
while lev<=depth
 
    %--------------------------------------------------------------
    % UNIFORM (UNIFORM QUADTREE)
    numbox = box_numbers(lev+1,1:numlev(lev+1));
    centB = cent(numbox,:);
    cent_intB = cent_int(numbox,:);         
            
    for i=1:numlev(lev+1)
        iB = cent_intB(i,1); jB = cent_intB(i,2); 
        
        % We use the sparse matrix Sindex{lev+1} to create an inverse that
        % maps (iB,jB) -> box number
        Sindex{lev+1}(iB,jB) = numbox(i);        
    end
            
    if lev<depth
        [Sindex,cent,cent_int,C,P,box_numbers,num_bisec,numnodes,numlev,levbox]=bisect_level(centB,numbox,levbox,lev,Sindex,cent,cent_int,C,P,box_numbers,num_bisec,numnodes,numlev);
    else
        % If bisection is not necessary, B is a leaf, and as such, has
        % no children. 
        C(numbox,1:2) = zeros(numlev(depth+1),nchild);
    end    

    
    % The number of boxes at all levels higher than lev might have changed
    % and so it is recomputed by multiplying the number of bisections by
    % nchild
    numlev(2:lev+2) = nchild*num_bisec(1:lev+1);
    
    % If no boxes were bisected, the loop must terminate. Otherwise, we go
    % down one level. 
    if num_bisec(lev+1)==0
        maxlev=lev;
        lev = depth + 1;
    else
        lev = lev + 1;
    end

end

leaves = zeros(2^depth,2); 

for lev=depth
    for k=1:numlev(lev+1)
          
    % We now use the array box_numbers to compute each box's center
    % and integer coordinates (iB,jB)
    numbox = box_numbers(lev+1,k);
            
    % Add boxes to the leaves array
    if C(numbox,1) == 0
        num_leaf = num_leaf+1;
        leaves(num_leaf,:) = [numbox , lev];
    end 
    end
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% SKELETON INFORMATION 
    
% Array and Cell Initialization
I_src = cell(nboxes,1); 
I_sk    = I_src; 
n_skel = zeros(nboxes,1); 
nbox = size(cent,1); 
k = zeros(depth+1,1); 
J = cell(depth+1,1);
    
layers = params.layers; %Max number of bdry layers for box skeletons, usually between 1 and 4.     
nmx = params.max_particles; 
    
I_source = Build_I_source_unif(X_source(1:nmx:N,:),N,nbox,nmx,depth,Sindex{depth+1});

% Parameters
h = params.h;
kh = params.kh; 
    
%Upward pass in which we build I_box and arrays associated with
%skeleton information
for lev=depth:-1:0
   for i = 1:numlev(lev+1)
       Nbox = box_numbers(lev+1,i);
       cB = cent(Nbox,:); 
          
   if mod(lev,2) == 0
       scl = lev/2;
       colrad = 2^(-scl); 
       rowrad = colrad; 
   else
       scl = (lev-1)/2;
       rowrad = 2^(-scl); 
       colrad = rowrad/2;
   end
          
   % If the box is a leaf, we pick skeleton points from I_source.
   % Note that it is assumed that max_particles <= (2*layers-1)^2 in 
   % case of a grid (meaning leaf boxes only have that many layers).
          
   if C(Nbox,1) == 0
       I_src{Nbox} = I_source(Nbox,I_source(Nbox,:)>0); 
   else
       I_src{Nbox} = [I_sk{C(Nbox,1)} I_sk{C(Nbox,2)} ]; 
   end

   % Size of I_src
   n_skel(Nbox) = length(I_src{Nbox}); 
   X_src = X_source(I_src{Nbox},:); 
   
   % Build skeleton sets (boundary layers
   if i==1 
      fprintf('\n Skeletons at level %d \n',lev); 
   end
              
   % We take skeleton points from boundary layers                         
   if i==1
       [J{lev+1},k(lev+1)] = LOCAL_Preset_Skeleton_Info(X_src,cB,lev,colrad,rowrad,h,layers,n_skel(Nbox),kh);              
   end
                  
   I_sk{Nbox} = I_src{Nbox}(J{lev+1}(1:k(lev+1))); 
   
   end
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
    
% SAVE TREE DATA into the struct INFO_TREE
nboxes = size(cent,1); 

for NB = nboxes:-1:1
    TREE.BOX(NB) = struct('cent',cent(NB,:),'parent',P(NB),'child',C(NB,:),'levbox',levbox(NB),'I_src',I_src{NB},...
                    'n_skel',n_skel(NB),'I_sk',I_sk{NB});
end
            
for lev = maxlev:-1:0
    TREE.LEV(lev+1) = struct('k',k(lev+1),'J',J{lev+1}); 
end

TREE.numlev = numlev; 
TREE.leaves = leaves; 
TREE.box_numbers = box_numbers;
TREE.depth = depth; 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [I_source] = Build_I_source_unif(Xi,N,nbox,nmx,depth,S)
    ind = 1:nmx:N; 
    nX = size(Xi,1); 
    lev = depth; 
    if mod(lev,2) == 0
        scl = lev/2;
        % level is even, boxes are square, scl = lev/2;
        coord = ceil((2^(scl-1))*(Xi + ones(nX,2)));
        %correcting for points with coordinates == -1. 
        coord = coord + (coord==0); 
    else
         scl = (lev-1)/2;
         % level is odd, boxes are rectangles, scl = (lev-1)/2;
         coord(:,1) = ceil((2^(scl))*(Xi(:,1) + ones(nX,1)));
         coord(:,2) = ceil((2^(scl-1))*(Xi(:,2) + ones(nX,1)));
         %correcting for points with coordinates == -1. 
         coord = coord + (coord==0); 
    end
    
    I_source = zeros(nbox,nmx); 
    
    for i=1:nX
        Nbox = S(coord(i,1),coord(i,2)); 
        I_source(Nbox,:) = ind(i):ind(i)+nmx-1; 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [J,k] = LOCAL_Preset_Skeleton_Info(X_src,cB,lev,colrad,rowrad,h,layers,n_skel,kh)
global BOX 

 % BOX DIMENSIONS
 bdup1 = cB(1) + colrad - (layers)*h; 
 bddn1 = cB(1) - colrad + (layers)*h; 
 bdup2 = cB(2) + rowrad - (layers)*h; 
 bddn2 = cB(2) - rowrad + (layers)*h;
              
 % SKELETON POINTS - BOUNDARY LAYERS
 Idx = (1:n_skel)'; 
 sk_ind = (X_src(:,1) >= bdup1 | X_src(:,1) <= bddn1 | X_src(:,2) >= bdup2 | X_src(:,2) <= bddn2); 
              
 % THE REST OF THE POINTS ARE ON THE INTERFACE
 %res_ind = (1:n_skel);
 %res_ind(J(1:k)) = [];  
 
 if mod(lev,2)==0
     norm_int = 1; tan_int = 2; rad = rowrad; th0 = pi/2; 
 else
     norm_int = 2; tan_int = 1; rad = colrad; th0 = 0; 
 end
 
 lambda = 2*pi/(max(kh,h)); wvls = ceil(rad/(lambda));
 
 if kh>30 && wvls>=1
     % Points not on interface, from previous levels
    nint_ind = (~sk_ind)&(X_src(:,norm_int) > cB(norm_int) + (layers)*h | X_src(:,norm_int) < cB(norm_int) - (layers)*h ); 
    sk_or_nint = (sk_ind | nint_ind); 
        
    % Sort interfacial points
    X_rs = X_src(~sk_or_nint,:); 
    [~,J1] = sort(X_rs(:,norm_int)); 
    [~,J2] = sort(-X_rs(J1,tan_int));  
    Jr = J1(J2); 
       
    % Add points per # of wavelengths on interface
    L = length(Jr); st = floor(2*L/((wvls+2))); 
    I_rs = Idx(~sk_or_nint); I_rs = I_rs(Jr); 
    
    extra_ind = floor(st/2)+1:st:L - floor(st/2) + 1; 
    md = mod(extra_ind - 1,2*layers); 
    extra_ind = [extra_ind min((extra_ind + (2*layers - 1 - 2*md)),L)];
    
    %if pB == 0
    %   extra_ind = (L+1) - extra_ind;  
    %end
    
    sk_ind(I_rs(extra_ind)) = true;
    I_rs(extra_ind) = [];
    
    % Add extra points and re-sort
    I_nint = Idx(nint_ind);
    X_nint = X_src(I_nint,:); 
    X_rs = [X_src(I_rs,:) ; X_nint]; 
    I_rs = [I_rs ; I_nint];
    [~,J1] = sort(X_rs(:,norm_int)); 
    [~,J2] = sort(-X_rs(J1,tan_int));  
    Jr = J1(J2); 
    I_rs = I_rs(Jr); 
 else
    % Sort interfacial points
    X_rs = X_src(~sk_ind,:); 
    [~,J1] = sort(X_rs(:,norm_int)); 
    [~,J2] = sort(-X_rs(J1,tan_int));  
    Jr = J1(J2);   
    I_rs = Idx(~sk_ind); I_rs = I_rs(Jr); 
 end
 
 % Sort skeleton points
 I_sk = Idx(sk_ind); 
 X_sk = X_src(sk_ind,:); 
 [th,rho] = cart2pol((X_sk(:,1) - cB(1)),X_sk(:,2) - cB(2));
 th(th<th0) = th(th<th0) + 2*pi; 
 [~,J1] = sort(rho);
 [~,J2] = sort(th(J1));
 Js = J1(J2);
 I_sk = I_sk(Js); 
 
 % J = [I_sk I_rs]; 
 k = size(I_sk,1); 
 J = [I_sk ; I_rs];  

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [S,cent,cent_int,C,P,box_numbers,num_bisec,numnodes,numlev,levbox] = bisect_level(centB,numbox,levbox,lev,S,cent,cent_int,C,P,box_numbers,num_bisec,numnodes,numlev)
    % This function performs the bisection of box B (with center centB, at
    % level lev with coordinates (iB,jB)). 
    
    nchild = 2^(lev+1);
    if mod(lev,2) == 0
        % Vertical Bisection
        cta = [-0.5 0 ; 0.5 0];
        scl = lev/2;
    else
        % Horizontal Bisection
        cta = [0 0.5 ; 0 -0.5];
        scl = (lev-1)/2;
    end
    
    cta = repmat(cta,2^lev,1); 
    
    
    % Update to data structures (perform the bisection). We assign the last
    % 2 numbers to the new children boxes. 
    % Centers
    childidx = (1:nchild) + numnodes*ones(1,nchild);
    
    cent(childidx,:) = centB(ceil((1:nchild)./2),:)+cta*(1/2)^(scl);
    % Levels
    levbox(childidx) = (lev+1)*ones(1,nchild);
    % Now our box has nchild children
    C(numbox,1:2) = reshape(childidx,2,2^lev)';
    % These children are, as of now, leaves (no children of their own)
    C(childidx,:) = zeros(nchild,2); 
    % Box Numbers
    box_numbers(lev+2,1:nchild) = childidx;
    % nchild more on this level
    numlev(lev+2) = nchild;
    % Parent info for children
    P(childidx) = floor(0.5*childidx);
    
    % We add the children's integer coordinates to S(i,j)
    for ch=1:nchild
        if mod(lev+1,2) == 0
            iC = ceil((2^(scl))*(cent(numnodes+ch,1)+1));
            jC = ceil((2^(scl))*(cent(numnodes+ch,2)+1));
        else
            iC = (ceil((2^(scl+1))*(cent(numnodes+ch,1)+1))+1)/2;
            jC = (ceil((2^(scl))*(cent(numnodes+ch,2)+1))+1)/2;
        end
        
        cent_int(ch + numnodes,:) = [iC jC]; 
        S{lev+2}(iC,jC) = numnodes+ch;
    
    end
    
    % Finally, the total number of nodes increases by nchild, and the
    % number of bisections by one. 
    numnodes = numnodes + nchild;
    num_bisec(lev+1) = 2^lev;

end