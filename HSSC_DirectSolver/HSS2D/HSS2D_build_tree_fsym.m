function [INFO_TREE]=HSS2D_build_tree_fsym(params)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%     FUNCTION CALL:
%         [INFO_TREE] = HSS2D_build_tree_fsym(params)
%
%     DESCRIPTION:
%         This function does most of the work in compressing A as HSS2D. It constructs the
%         binary tree with the index sets and interpolation operators for each box. It is
%         assumed that the associated kernel is symmetric (fsym).
%
%         The HSS2D data structure is a binary tree (essentially a quadtree) subdivision of a
%         planar surface of maximum depth params.depth. It is constructed by recursively
%         alternating between horizontal and vertical bisection of boxes. If params.bisec_rule =
%         'maxpart', the tree is adaptive and level-restricted. This means that a box is refined
%         only if the number of source points >= params.max_particles and its level is
%         <= params.depth.
%
%         The points associated with a box I_src are taken to be the union of the skeleton
%         points of children. These points are split into skeleton I_sk and residual I_rs sets
%         in a way depending on params. The method described in the paper takes skeleton points
%         to be those along the boundary of the box and the residual points to be those along
%         the interface of the children (params.skel_rule = 'bdry').
%
%         The interpolation matrices T that comprise the left and right factors L and R in the
%         telescoping factorization are stored for each box. The fast method of computing them
%         is described in Algorithm 3 and Remark 3.4 of the paper. This corresponds to using
%         preset boundary layer skeletons (params.skel_rule = 'bdry'), and computing T as
%         lowrank (params.INTERform = 'lowrank'). In this case, an HSS1D least squares and
%         randomized ID speedup are used for any box with more points than params.k_cut (default
%         2500).
%
%         Because special care is needed for 2D compression, the computation of the self (leaf)
%         and sibling interaction (non-leaf) matrices that comprise D in the telescoping
%         factorization is held off until later. They are implicitly stored by the index sets.
%
%
%     INPUT:
%         params (see HSS_tree_parameters.m) HSS2D parameters.
%
%
%     OUTPUT:
%         INFO_TREE       <struct>        The data structure for HSS2D. It stores the
%                                             interpolation matrices and index sets in a binary
%                                             tree.
%         .depth          <int>           Maximum level of the tree.
%         .numlev         <depthx1 int>   Number of boxes per level.
%         .box_numbers    <sparse int>    Array that stores the box numbers by level in order to
%                                             traverse the tree.
%         .leaves         <mx2 int>       Index and level of each leaf box.
%
%         .BOX(i) (info for box, including index sets and ID matrices)
%             .cent   <1x2 float>         Box center coordinates.
%             .levbox <int>               Box level in the tree.
%             .parent <int>               Box parent index (0 for root node).
%             .child  <1x2 int>           Box's children indices (non leaf boxes).
%             .I_src  <n_skelx1 int>      All source points in the box.
%             .I_sk   <kx1 int>           Skeleton points in the box sorted by polar angle
%                                             around the box center.
%             .n_skel <int>               Size of I_src.
%             .k      <int>               Size of skeleton I_sk.
%             .J      <n_skelx1 int>      [I_sk , I_rs] Index with upwards skeleton first,
%                                             residual points second. Residual points are sorted
%                                             along the children interface.
%             .T      <kxn_skel-k float> Interpolation matrix that is stored either densely or
%                                             as a low-rank (T.U,T.V). MATLAB functions
%                                             interpolation_operator_*.m are used for preset
%                                             skeletons. Otherwise ID.m is used.
%
%         .LEV(i) (info for level i, for TI kernels)
%             If params.transinv == True, J and T are stored here rather than in .BOX()
%




depth = params.depth; 
X_source = params.X_source;
N = size(X_source,1);

% Array Initialization
nboxes = 2^(depth+1)-1;
cent = zeros(nboxes,2); 
cent_int = cent; 
C = cent; P = cent(:,1); 
levbox = P;
count = P; 
I_src = cell(nboxes,1); 
neighbors = I_src; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
count(1) = N; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (I) BINARY TREE CONSTRUCTION

% We traverse the tree while the level of the boxes considered is not
% greater than the maximum assigned depth, and as long as bisection is
% required. 
% When required, a bisection function is used to ensure that the tree
% is bisected adaptively and/or is level restricted. 
% On the first pass, the necessary arrays are created for each box going 
% down the binary tree. 

lev=0;
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
    
    scl = floor(lev/2); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Integer Coordinates for Maximum Particle (Adaptive) build
    if strcmp(params.bisec_rule,'maxpart')
        if lev>0
        % The array coord gives us integer coordinates (i,j) (i,j in
        % {1,..,2^lev}) for each point in X_source, and we record the number of
        % particles in each box in count. 
        if mod(lev,2) == 0
            % level is even, boxes are square, scl = lev/2;
            tmp = (2^(scl-1))*(X_source+1);
            coord = ceil(tmp);
        else
            % level is odd, boxes are rectangles, scl = (lev-1)/2;
            tmp(:,1) = (2^(scl))*(X_source(:,1)+1);
            tmp(:,2) = (2^(scl-1))*(X_source(:,2)+1);
            coord = ceil(tmp);
        end
        
        %correcting for points with coordinates == -1. 
        coord = coord + (coord==0);
         
        n2 = 2^scl; 
        tag     = coord(:,1)+n2*coord(:,2)-n2; 
        numbox  = box_numbers(lev+1,1:numlev(lev+1));
        tag_box = cent_int(numbox,1) + n2*cent_int(numbox,2)-n2; 
        ind = repmat(numbox(tag),numlev(lev+1),1)==repmat(box(tag_box).',1,N);
        Idx = (1:N).';
               
        % Adding points on each box
        count(numbox) = sum(ind.').'; 
        for i=1:numlev(lev+1)
            I_src{numbox(i)} = Idx(ind(i,:)); 
        end
        
        Sindex{lev+1}(tag_box) = numbox; 
        
        else
           count(1) = N; 
           I_src{1} = (1:N).';
        end
    end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We consider each box on the current level, and we bisect using the function bisect_node 
 
% MAXPART (MAXIMUM PARTICLES CRITERION)
        
    if strcmp(params.bisec_rule,'maxpart')
        for k = 1:numlev(lev+1)
            % Box information
            numbox = box_numbers(lev+1,k);
            centB = cent(numbox,:);    
        
            %------------------------------------------------------------------
            % if the number of particles is greater than n_max, we bisext into 
            % 2 new child boxes and add information to our arrays as needed
            if count(numbox) > params.max_particles && lev<depth
            
            [cent,cent_int,C,P,box_numbers,num_bisec,numnodes,numlev,levbox]...
                =bisect_node(centB,numbox,levbox,lev,cent,cent_int,C,P,box_numbers,num_bisec,numnodes,numlev);
            else
            % If bisection is not necessary, B is a leaf, and as such, has
            % no children. 
            C(numbox,1:nchild) = zeros(1,nchild);
            end
        end
    else
        %--------------------------------------------------------------
        % UNIFORM (UNIFORM BINARY TREE)
        if strcmp(params.bisec_rule,'uniform')
            numbox = box_numbers(lev+1,1:numlev(lev+1));
            centB = cent(numbox,:);
            cent_intB = cent_int(numbox,:);  
            
            n2 = size(Sindex{lev+1},1);
            iB = cent_intB(:,1); jB = cent_intB(:,2); 
            ind = iB + n2*jB - n2;
            Sindex{lev+1}(ind) = numbox;   
            
            if lev>1
            
            shift = [-1 -1;-1 0;-1 1;0 -1;0 1;1 -1;1 0;1 1]; 
            
            for i=1:numlev(lev+1)
                MI = max(iB); MJ = max(jB);
                % Integer coordinates of neighbor box centers
                crd_ngh = repmat(cent_int(numbox(i),:),8,1) + shift; 
                valid1 = (crd_ngh(:,1)>0 & crd_ngh(:,1)<=MI); 
                valid2 = (crd_ngh(:,2)>0 & crd_ngh(:,2)<=MJ);
                
                if params.dim ==3 & strcmp(params.type,'torus')
                    crd_ngh(~valid1,1) = mod(crd_ngh(~valid1,1),MI);   
                    crd_ngh(~valid2,2) = mod(crd_ngh(~valid2,2),MJ);    
                    crd_ngh(crd_ngh(:,1)==0,1) = MI; 
                    crd_ngh(crd_ngh(:,2)==0,2) = MJ; 
                else
                    crd_ngh = crd_ngh(valid1&valid2,:);
                end
                
                % Obtain neighbor box numbers using Sindex
                neighbors{numbox(i)} = full(Sindex{lev+1}(crd_ngh(:,1)+n2*crd_ngh(:,2)-n2));   
                neighbors{numbox(i)} = unique(neighbors{numbox(i)});   
            end
            else
               if lev>0
                   neighbors{2} = 3;
                   neighbors{3} = 2;   
               else
                    neighbors{1} = [];    
               end
            end
            
            if lev<depth
                [cent,cent_int,C,P,box_numbers,num_bisec,numnodes,numlev,levbox]...
                    =bisect_level(centB,numbox,levbox,lev,cent,cent_int,C,P,box_numbers,num_bisec,numnodes,numlev);
            else
                % If bisection is not necessary, B is a leaf, and as such, has
                % no children. 
                C(numbox,1:2) = zeros(numlev(depth+1),nchild);
            end    
            
        end
    end
    
    % If no boxes were bisected, the loop must terminate. Otherwise, we go
    % down one level. 
    if num_bisec(lev+1)==0
        maxlev=lev;
        lev = depth + 1;
    else
        lev = lev + 1;
    end 
end

nboxes = size(cent,1); 
Id_box = 1:nboxes; 
leaf_ind = C(:,1)==0 & C(:,2)==0;
leaves(:,1) = Id_box(leaf_ind); 
leaves(:,2) = levbox(leaf_ind); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SKELETON INFORMATION (if DOskeletons is 1)
if (params.DOskeletons)
    
    % Parameters 
    layers = params.layers;
    nmx = params.max_particles; 
    h = params.h; 
    knl_par = params; knl_par.X_source = []; 
    q = zeros(depth+1,1); 
    width = layers*h; 
    
    % Array and Cell Initialization
    if strcmp(params.bisec_rule,'uniform')
        I_src = cell(nboxes,1); 
        I_sk    = I_src; 
        n_skel = zeros(nboxes,1); 
        Xp = I_src;  
        if params.transinv == 0
            k = n_skel; T = I_sk; J = T;
        else
            k = zeros(depth+1,1); T = cell(depth+1,1); J = T;
        end
        
        % Build I_src{Nbox}. For the case of an adaptive tree, it should be
        % done per level in leaves. 
        [I_src,~] = Build_I_source_unif(X_source(1:nmx:N,:),N,nboxes,nmx,depth,Sindex{depth+1});
    end
    
    %Upward pass in which we build I_box and arrays associated with
    %skeleton information
    for lev=maxlev:-1:0
       scl = floor(lev/2);
       rowrad = 2^-scl; 
       mdl = mod(lev,2); 
       colrad = rowrad/(2^mdl); 
       rad = min(rowrad,colrad); 
       
       if mdl == 0  
           th0 = pi/2;   
           th1 = -pi/2;  
       else
           th0 = 0; 
           th1 = -pi; 
       end
       
       norm_int = 1+mdl; tan_int = 2-mdl;
       
       % Layers of points on box proxies
       if lev>0 | params.m<1       
          if strcmp(params.flag_pot,'SL_L_2D') 
            lay_prox = layers+1; 
          else
            lay_prox = 2*round(rad/h);   
          end
       else
          lay_prox = layers;     
       end   
     
       fprintf('\n layers on proxy = %d',lay_prox)
        
       for i = 1:numlev(lev+1)
          Nbox = box_numbers(lev+1,i);
          cB = cent(Nbox,:); 
          
          
          % If the box is a leaf, we pick skeleton points from I_source.
          % If it's not a leaf, we merge skeleton sets from children. 
          if C(Nbox,1) > 0 & C(Nbox,2)>0
              I_src{Nbox} = [I_sk{C(Nbox,1)} I_sk{C(Nbox,2)} ];
          elseif C(Nbox,1)>0 | C(Nbox,2)>0
              child = max(C(Nbox,1),C(Nbox,2)); 
              I_src{Nbox} = I_sk{child}; 
          end
          
          % Number of points inside box 
          n_skel(Nbox) = length(I_src{Nbox}); 
          X_src = X_source(I_src{Nbox},:); 
          
          % Skeleton set parameters: 
          % Preset skeletons     - boundary layers or via a dyadic cube
          % decomposition, with the option of adding extra points with a
          % randomized ID. 
          % Non Preset skeletons - determined via an ID on K[X_proxy,X_sk] 
          if i==1
              knl_par.proxy_rule = params.proxy_rule; 
            if lev<depth-1 || (params.kh == 0)                
                
                if strcmp(params.skel_rule,'ID')
                    params.preset = 0;
                    add_extra_pts = 0;
                else
                    params.preset = 1; 
                    add_extra_pts = params.skel_extra_pts;     
                end
            else
                params.preset = 0;
                 add_extra_pts = 0; 
            end
            
            if lev==0 & params.m==1   
                add_extra_pts = 0;     
            end
              
            fprintf('\n Skeletons at level %d, preset = %d, add_extra_pts = %d',lev,params.preset,add_extra_pts); 
          end    
              
          if params.preset == 0
              % NON PRESET Skeletons (ID) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                if params.transinv == 0
                    % Non Translation Invariant (T,J,k)(Nbox)
                        
                    % Build proxy points
                    if strcmp(params.proxy_rule,'neigh')
                        tr_par.X_source = X_source; 
                        tr_par.neighbors = neighbors{Nbox}; 
                        tr_par.C = C;
                        tr_par.I_sk = I_sk; 
                        tr_par.I_src = I_src; 
                        
                        if lev>0 & lev<depth-1
                            X_ext = LOCAL_Build_Proxy(cB,colrad,rowrad,lev,n_skel(Nbox),lay_prox,knl_par,tr_par); 
                            size(X_ext)      
                        else
                            knl_par.proxy_rule = 'bdry';   
                            if lev>=0
                                lay = lay_prox;    
                                X_ext = LOCAL_Build_Proxy(cB,colrad,rowrad,lev,n_skel(Nbox),lay,knl_par); 
                            else
                                lay = layers;   
                                X_ext = LOCAL_Build_Proxy(cB,1-width,1-width,lev,n_skel(Nbox),lay,knl_par);
                            end      
                        end
                    
                    else
                        
                        X_ext = LOCAL_Build_Proxy(cB,colrad,rowrad,lev,n_skel(Nbox),lay_prox,knl_par);
                    %{
                    if lev<depth-1 && lev>0      
                        Inside = (abs(X_ext(:,1))<1 & abs(X_ext(:,2))<1);    
                        X_ext = X_ext(Inside,:); 
                    else
                        Inside = (abs(X_ext(:,1))<1+width & abs(X_ext(:,2))<1+width);    
                        X_ext = X_ext(Inside,:);    
                    end
                    %}
                    end     
              
                    % Skeleton Info using the ID
                    [T{Nbox},J{Nbox},k(Nbox)] = LOCAL_Non_Preset_Skeleton_Info(X_ext,X_src,cB,lev,colrad,rowrad,knl_par);
                        
                    I_sk{Nbox} = I_src{Nbox}(J{Nbox}(1:k(Nbox)));    
                    
                    Xp{Nbox} = X_ext; 
                else
                    % Translation Invariant (T,J,k)(level)
                    if i==1
                        X_ext = LOCAL_Build_Proxy(cB,colrad,rowrad,lev,n_skel(Nbox),lay_prox,knl_par);
              
                        % Skeleton Info using the ID
                        [T{lev+1},J{lev+1},k(lev+1)] = LOCAL_Non_Preset_Skeleton_Info(X_ext,X_src,cB,lev,colrad,rowrad,knl_par);   
                    end
               
                I_sk{Nbox} = I_src{Nbox}(J{lev+1}(1:k(lev+1))); 
                end
          else
              % PRESET Skeletons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
              if params.transinv == 0 
                 if mod(Nbox,2)==0
                     th = th0; 
                 else
                     th = th1; 
                 end
                  
                 % Non Translation Invariant (T,J,k)(Nbox)
                [J{Nbox},k(Nbox)] = LOCAL_Preset_Skeleton_Info(X_src,cB,lev,colrad,rowrad,knl_par);
                
                % Build proxy points  
                if strcmp(params.proxy_rule,'neigh')
                    if lev>0
                    NN = length(neighbors{Nbox}); 
                    IN = [];
                    for j=1:NN
                        ngh = neighbors{Nbox}(j);    
                        if C(ngh,1) > 0 & C(ngh,2)>0
                            IN = [IN I_sk{C(ngh,1)} I_sk{C(ngh,2)} ];
                        elseif C(ngh,1)>0 | C(ngh,2)>0
                            child = max(C(ngh,1),C(ngh,2)); 
                            IN = [IN I_sk{child}];
                        else
                            IN = [IN I_src{ngh}]; 
                        end
                    end
                    
                    X_ext = X_source(IN,:); 
                    %D_ext = abs(min(rad-abs(X_ext(:,1)-cB(1)),rad-abs(X_ext(:,2)-cB(2)))); 
                    %X_ext = X_ext(D_ext<2*lay_prox*h,:);                  
                    
                    if params.dim==3   
                       isk = J{Nbox}(1:k(Nbox));
                       X_ext = shape_3D(X_ext,params,params.type);  
                       CB    = shape_3D(cB,params,params.type); 
                       Z_sk  = shape_3D(X_src(isk,:),params,params.type); 
                       brad(1) = max(abs(Z_sk(:,1)-CB(1)));    
                       brad(2) = max(abs(Z_sk(:,2)-CB(2)));
                       brad(3) = max(abs(Z_sk(:,3)-CB(3)));
                       D_ext = abs(min([brad(1)-abs(X_ext(:,1)-CB(1)),...
                                        brad(2)-abs(X_ext(:,2)-CB(2)),...
                                        brad(3)-abs(X_ext(:,3)-CB(3))].')); 
                       br = min(brad); 
                       %[~,Js] = sort(D_ext);   
                       ps = min(n_skel(Nbox),size(X_ext,1)); 
                       ind_ex = D_ext<2*br;   
                       while(sum(ind_ex)<ps)    
                           br = 2*br;   
                           ind_ex = D_ext<1.5*br;   
                           sum(ind_ex)
                           display(br)
                       end
                       X_ext = X_ext(ind_ex,:);      
                    end
                    
                    size(X_ext)   
                    
                    else
                        knl_par.proxy_rule = 'bdry';   
                        X_ext = LOCAL_Build_Proxy(cB,colrad,rowrad,lev,n_skel(Nbox),layers,knl_par);   
                    end
                    
                else
                    if lev==-1   
                        lay = layers; 
                        knl_par.proxy_rule = 'bdry';   
                        X_ext = LOCAL_Build_Proxy(cB,1-width,1-width,lev,n_skel(Nbox),lay,knl_par);
                    else
                        X_ext = LOCAL_Build_Proxy(cB,colrad,rowrad,lev,n_skel(Nbox),lay_prox,knl_par);
                    end   
                    %{
                    if lev>0  
                        Inside = (abs(X_ext(:,1))<1 & abs(X_ext(:,2))<1);    
                        X_ext = X_ext(Inside,:); 
                    else
                        Inside = (abs(X_ext(:,1))<1+width & abs(X_ext(:,2))<1+width);    
                        X_ext = X_ext(Inside,:);     
                    end
                        %}
                end  
                
                Xp{Nbox} = X_ext; 
              
                % Compute Interpolation Operator K[X_ext,X_sk]\K[X_ext,X_rs]
                T{Nbox} = LOCAL_Interpolation_Operators(cB,lev,maxlev,params,X_ext,X_src,J{Nbox},k(Nbox),q);
                
                % If we require more points (and k<n), an update of {T,J,k} 
                %is performed using a randomized ID. 
                if add_extra_pts & (n_skel(Nbox) - k(Nbox)) > 0
                    [T{Nbox},J{Nbox},k(Nbox)] = LOCAL_ID_rand_update(T{Nbox},J{Nbox},k(Nbox),X_ext,X_src,knl_par); 
                    
                    % Orient skeletons after adding points
                    J_sk = J{Nbox}(1:k(Nbox));
                    J_rs = J{Nbox}(k(Nbox)+1:end); 
                
                    Js = LOCAL_orient_skeleton(X_src(J_sk,:),cB,th);    
                    %Jr = LOCAL_orient_residual(X_src(J_rs,:),norm_int,tan_int);
                    Jr = 1:length(J_rs);    
                
                    J{Nbox} = [J_sk(Js) ; J_rs(Jr)];   
                    
                    if isempty(T{Nbox}.V)
                       T{Nbox}.U = T{Nbox}.U(Js,Jr);  
                    else
                       T{Nbox}.U = T{Nbox}.U(Js,:); 
                       T{Nbox}.V = T{Nbox}.V(Jr,:); 
                    end    
                end
                
                % Define I_sk (global indices) 
                I_sk{Nbox} = I_src{Nbox}(J{Nbox}(1:k(Nbox))); 
              else
                % Translation Invariant (T,J,k)(level)  
                if i==1
                    [J{lev+1},k(lev+1)] = LOCAL_Preset_Skeleton_Info(X_src,cB,lev,colrad,rowrad,knl_par);  
                    
                    % Build proxy points 
                    X_ext = LOCAL_Build_Proxy(cB,colrad,rowrad,lev,n_skel(Nbox),lay_prox,knl_par);
                    
                    % Compute Interpolation Operator K[X_ext,X_sk]\K[X_ext,X_rs]
                    
                    T{lev+1} = LOCAL_Interpolation_Operators(cB,lev,maxlev,params,X_ext,X_src,J{lev+1},k(lev+1),q);
                    
                    % If we require more points (and k<n), an update of {T,J,k} 
                    %is performed using a randomized ID. 
                    if add_extra_pts & (n_skel(Nbox) - k(lev+1)) > 0
                        [T{lev+1},J{lev+1},k(lev+1)] = LOCAL_ID_rand_update(T{lev+1},J{lev+1},k(lev+1),X_ext,X_src,knl_par); 
                        
                        % Orient skeletons after adding points
                        J_sk = J{lev+1}(1:k(lev+1));
                        J_rs = J{lev+1}(k(lev+1)+1:end); 
                
                        Js = LOCAL_orient_skeleton(X_src(J_sk,:),cB,th0); 
                        Jr = LOCAL_orient_residual(X_src(J_rs,:),norm_int,tan_int);  
                
                        J{lev+1} = [J_sk(Js) ; J_rs(Jr)];   
                    
                        if isempty(T{lev+1}.V)
                            T{lev+1}.U = T{lev+1}.U(Js,Jr);  
                        else
                            T{lev+1}.U = T{lev+1}.U(Js,:); 
                            T{lev+1}.V = T{lev+1}.V(Jr,:); 
                        end
                    end
                   
                end
                  
                I_sk{Nbox} = I_src{Nbox}(J{lev+1}(1:k(lev+1))); 
              end
          end
         
         % Rank of interpolation operators 
         if params.transinv==0 && lev>0 && size(T{Nbox},2)>0
            if i==1 
                q(lev+1) = size(T{Nbox}.U,2);  
            else
                q(lev+1) = max(q(lev+1),size(T{Nbox}.U,2)); 
            end
         end
          
       end
       
       if params.transinv==1 && lev>0   
          q(lev+1) = size(T{lev+1}.U,2);  
       end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% SAVE TREE DATA into the struct INFO_TREE
nboxes = size(cent,1); 

if strcmp(params.bisec_rule,'maxpart')
    if params.DOskeletons
        for NB = nboxes:-1:1
            if C(NB,1) == 0
                INFO_TREE.BOX(NB) = struct('cent',cent(NB,:),'parent',P(NB),'child',C(NB,:),'levbox',levbox(NB),'I_src',I_src{NB},...
                'count',count(NB),'n_skel',n_skel(NB),'k',k(NB),'J',J{NB},'T',T{NB},'I_sk',...
                I_sk{NB});
            else
                INFO_TREE.BOX(NB) = struct('cent',cent(NB,:),'parent',P(NB),'child',C(NB,:),'levbox',levbox(NB),'I_src',I_src{NB},...
                'count',count(NB),'n_skel',n_skel(NB),'k',k(NB),'J',J{NB},'T',T{NB},'I_sk',...
                I_sk{NB});
            end
        end
    else
        for NB = nboxes:-1:1
            INFO_TREE.BOX(NB) = struct('cent',cent(NB,:),'parent',P(NB),'child',C(NB,:),'levbox',levbox(NB));  
        end
    end
    INFO_TREE.numlev = numlev; 
    INFO_TREE.leaves = leaves; 
    INFO_TREE.Sindex = Sindex;
    INFO_TREE.box_numbers = box_numbers;
    INFO_TREE.depth = maxlev; 
else
    if params.DOskeletons
        if params.transinv            
            for NB = nboxes:-1:1
                INFO_TREE.BOX(NB) = struct('cent',cent(NB,:),'parent',P(NB),'child',C(NB,:),'neighbors',neighbors{NB},'levbox',levbox(NB),'I_src',I_src{NB},...
                    'n_skel',n_skel(NB),'I_sk',I_sk{NB});
            end
            
            for lev = maxlev:-1:0
                INFO_TREE.LEV(lev+1) = struct('k',k(lev+1),'J',J{lev+1},'T',T{lev+1}); 
            end
            
        else
            for NB = nboxes:-1:1
                INFO_TREE.BOX(NB) = struct('cent',cent(NB,:),'parent',P(NB),'child',C(NB,:),'neighbors',neighbors{NB},'levbox',levbox(NB),'I_src',I_src{NB},...
                    'n_skel',n_skel(NB),'k',k(NB),'J',J{NB},'T',T{NB},'I_sk',I_sk{NB});       
            end
        end
    else
        for NB = size(cent,1):-1:1  
            INFO_TREE.BOX(NB) = struct('cent',cent(NB,:),'parent',P(NB),'child',C(NB,:),'neighbors',neighbors{NB},'levbox',levbox(NB));  
        end
    end
    INFO_TREE.numlev = numlev; 
    INFO_TREE.leaves = leaves; 
    INFO_TREE.box_numbers = box_numbers;
    INFO_TREE.depth = maxlev; 
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,J,k] = LOCAL_Non_Preset_Skeleton_Info(X_ext,X_src,cB,lev,colrad,rowrad,par)
% SUBSELECT POINTS FROM LEAF NODE MATRICES 

acc = par.acc; 

if par.dim == 3
    Z_src = shape_3D(X_src,par,par.type);    
    A21 = Kernel_Eval(X_ext,Z_src,par); 
else
    A21 = Kernel_Eval(X_ext,X_src,par); 
end 

[Ttmp, J ] = ID(A21, acc);
k = size(Ttmp,1);

if mod(lev,2)==0
    norm_int = 1; tan_int = 2; th0 = pi/2; 
else
    norm_int = 2; tan_int = 1; th0 = 0; 
end

I_sk = J(1:k); I_rs = J(k+1:end);
X_sk = X_src(I_sk,:); X_rs = X_src(I_rs,:);

%Sort residual points along interface
[~,J1] = sort(X_rs(:,norm_int)); 
[~,J2] = sort(-X_rs(J1,tan_int));  
Jr = J1(J2);   
I_rs = I_rs(Jr); 

%Sort skeleton points along a curve around cB
[th,rho] = cart2pol((X_sk(:,1) - cB(1)),X_sk(:,2) - cB(2));
th(th<th0) = th(th<th0) + 2*pi; 
[~,J1] = sort(rho);
[~,J2] = sort(th(J1));
Js = J1(J2);
I_sk = I_sk(Js);    

J = [I_sk I_rs];   
              
T.U = Ttmp(Js,Jr);  
T.V = eye(size(Ttmp,2));  

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [J,k] = LOCAL_Preset_Skeleton_Info(X_src,cB,lev,colrad,rowrad,params)
global BOX 
          
h = params.h; layers = params.layers; Ng = params.Ng; dim = params.dim; 
n_skel = size(X_src,1); 
Idx = (1:n_skel)';


if strcmp(params.skel_rule,'wdec')
    % Whitney decomposition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Skeleton points are defined on a dyadic Whitney decomposition of the
    % box into cubes well-separated from the box boundary.
    
    if dim==2 | dim==3
        % Whitney decomposition on a subset of the plane
        [~,sk_ind] = Whitney_decomposition(X_src,cB,colrad,rowrad,layers,h,Ng); 
    else
        % Whitney decomposition of a surface in 3D using an adaptive tree
        % and the intersection of the dyadic cubes with the surface patch. 
        
        c3 = (max(X_src(:,3))+min(X_src(:,3)))/2;
        r3 = (max(X_src(:,3))-min(X_src(:,3)))/2 +h;
        rad = [colrad rowrad r3];
        cent = [cB c3];    
    
        [~,sk_ind] = Whitney_decomposition_2Dsurf(X_src,cent,rad,layers,h);
    end
    
    if mod(lev,2)==0
        norm_int = 1; tan_int = 2; th0 = pi/2; 
    else
        norm_int = 2; tan_int = 1; th0 = 0; 
    end
    
    % Sort residual points along box interface
    X_rs = X_src(~sk_ind,:); 
    Jr = LOCAL_orient_residual(X_rs,norm_int,tan_int);
    I_rs = Idx(~sk_ind); 
    I_rs = I_rs(Jr); 
    
    I_sk = Idx(sk_ind); 
    X_sk = X_src(sk_ind,:); 
       
elseif strcmp(params.skel_rule,'rand')
    if mod(lev,2) == 0
       k0 = ceil((2/3)*n_skel); norm_int = 1; tan_int = 2; th0 = pi/2; 
    else
       k0 = ceil((3/4)*n_skel); norm_int = 2; tan_int = 1; th0 = 0;  
    end
    
    lambda = params.lambda;        
    
    k = 0; ks = k0;  
    while k0>1.1*k
        [X_sk,~,B] = LOCAL_random_sampler(X_src,ks,cB,colrad,rowrad,lambda,layers,params);
        k = size(X_sk,1); 
        if k0>1.1*k 
            ks = 2*ks;  
            lambda = 0.9*lambda; 
        end
    end
    
    sk_ind = B;     
    
    % Sort residual points along box interface
    X_rs = X_src(~sk_ind,:); 
    Jr = LOCAL_orient_residual(X_rs,norm_int,tan_int);
    I_rs = Idx(~sk_ind); 
    I_rs = I_rs(Jr);   
    
    I_sk = Idx(sk_ind); 
else
    % Boundary Layers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Preset initial approximation to skeleton: full outer boundary layers
    width = layers*h;
    %{
    if lev==1
        bdry = (abs(X_src(:,1) - cB(1)) > (colrad - width));
    else
    %}
        bdry = (abs(X_src(:,1) - cB(1)) > (colrad - width)) | (abs(X_src(:,2) - cB(2)) > (rowrad - width));
    %end
 
    % Optional extra "sparse" inner boundary layers
    sp_layers = 0;   
    extra_ind = false(size(bdry)); 
    for i=1:sp_layers
        %{
        if lev==1
            extra_ind_sp = (~bdry | ~extra_ind) & (abs(X_src(:,1) - cB(1)) > colrad - (width + i*h));      
        else
            extra_ind_sp = (~bdry | ~extra_ind) & ((abs(X_src(:,1) - cB(1)) > colrad - (width + i*h)) ...
                | (abs(X_src(:,2) - cB(2)) > rowrad - (width + i*h)));
        end
        %}
        
        itrue = Idx(extra_ind_sp ==  true);
 
        for j=1:i
            extra_ind_sp(itrue(j:sp_layers+1:end)) = false;
        end
    
        extra_ind = extra_ind | extra_ind_sp;
    end
    
    % Skeleton Set (boolean index)
    sk_ind = bdry | extra_ind; 
 
    if mod(lev,2)==0
         norm_int = 1; tan_int = 2; th0 = pi/2; 
    else
         norm_int = 2; tan_int = 1; th0 = 0; 
    end
    
    % Sort residual points along box interface
    X_rs = X_src(~sk_ind,:); 
    Jr = LOCAL_orient_residual(X_rs,norm_int,tan_int);
    I_rs = Idx(~sk_ind); I_rs = I_rs(Jr); 
    
    I_sk = Idx(sk_ind); 
    X_sk = X_src(sk_ind,:); 
end
 
% Sort skeleton points in a curve around center cB. Currently this is done
% by sorting in polar coordinates. 
Js = LOCAL_orient_skeleton(X_sk,cB,th0); 
I_sk = I_sk(Js); 

% Define k and J
k = size(I_sk,1); 
J = [I_sk ; I_rs];  

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J_sk = LOCAL_orient_skeleton(X_sk,cB,th0)

[th,rho] = cart2pol((X_sk(:,1) - cB(1)),X_sk(:,2) - cB(2));
th(th<th0) = th(th<th0) + 2*pi; 
[~,J1] = sort(rho);
[~,J2] = sort(th(J1));
J_sk = J1(J2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J_rs = LOCAL_orient_residual(X_rs,nrm,tan)

[~,J1] = sort(X_rs(:,nrm)); 
[~,J2] = sort(-X_rs(J1,tan));  
J_rs = J1(J2);   

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X_ext = LOCAL_Build_Proxy(cB,colrad,rowrad,lev,n_skel,lay_prox,params,tr_params)

layers = params.layers; h = params.h; 
if nargin<8
    tr_params = [];
end

if strcmp(params.proxy_rule,'wdec') 
    [X_ext,X_ext2] = build_box_proxy(cB,colrad,rowrad,lay_prox,h);    
    X_ext = [X_ext ; X_ext2];
                    
    [~,B] = Whitney_decomposition_1Dcurve(X_ext,[0 0],[1 1],cB,[colrad rowrad],2*layers,h);
    X_ext = X_ext(B,:); 
elseif strcmp(params.proxy_rule,'rand')
    [X_ext,X_ext2] = build_box_proxy(cB,colrad,rowrad,lay_prox,h);  
    X_ext = [X_ext ; X_ext2];
                                      
    if mod(lev,2) == 0
        k0 = ceil((2/3)*n_skel); 
    else
        k0 = ceil((3/4)*n_skel); 
    end
    
    k0 = min(k0,size(X_ext,1));  
    
    lambda = params.lambda;  
    k = 0; ks = k0; c = 0;   
    while k0>k & c<10  
        [X_ext,~,~] = LOCAL_random_sampler(X_ext,ks,cB,colrad,rowrad,lambda,layers,params);
        k = size(X_ext,1);   
        if k0>k 
            ks = 2*ks; 
            lambda = 0.9*lambda; 
        end
        c = c+1; 
    end

elseif strcmp(params.proxy_rule,'neigh')   
    C    = tr_params.C; 
    I_sk = tr_params.I_sk; 
    I_src = tr_params.I_src; 
    
    NN = length(tr_params.neighbors); 
    IN = [];
    for j=1:NN
        ngh = tr_params.neighbors{Nbox}(j);    
        if C(ngh,1) > 0 & C(ngh,2)>0
            IN = [IN I_sk{C(ngh,1)} I_sk{C(ngh,2)} ];
        elseif C(ngh,1)>0 | C(ngh,2)>0
            child = max(C(ngh,1),C(ngh,2)); 
            IN = [IN I_sk{child}];
        else
            IN = [IN I_src{ngh}]; 
        end
    end   
                    
    X_ext = tr_params.X_source(IN,:); 
    %D_ext = abs(min(rad-abs(X_ext(:,1)-cB(1)),rad-abs(X_ext(:,2)-cB(2)))); 
    %X_ext = X_ext(D_ext<lay_prox*h,:);       
                    
    if params.dim==3   
        X_ext = shape_3D(X_ext,params,params.type);  
        CB    = shape_3D(cB,params,params.type); 
        Z_src  = shape_3D(X_src,params,params.type); 
        brad(1) = max(abs(Z_src(:,1)-CB(1)));   
        brad(2) = max(abs(Z_src(:,2)-CB(2)));
        brad(3) = max(abs(Z_src(:,3)-CB(3)));
        D_ext = abs(min([brad(1)-abs(X_ext(:,1)-CB(1)) ...
                         brad(2)-abs(X_ext(:,2)-CB(2)) ...
                         brad(3)-abs(X_ext(:,3)-CB(3))].')); 
        br = min(brad); 
        %[~,Js] = sort(D_ext);   
        %ps = min(2*n_skel(Nbox),size(X_ext,1));             
        X_ext = X_ext(D_ext<1.5*br,:);            
    end     
else
                    
    [X_ext,X_ext2] = build_box_proxy(cB,colrad,rowrad,lay_prox,h);
    X_ext = [X_ext ; X_ext2]; 
end  
%{
if lev<=1                   
    Inside = abs(X_ext(:,2)-cB(2))<1;   
    X_ext = X_ext(Inside,:); 
end
%}  
                
if params.dim==3 
    X_ext = shape_3D(X_ext,params,params.type);   
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xr,Js,B] = LOCAL_random_sampler(X,k,cB,colrad,rowrad,lambda,lay,par)
rng('default');     % For older versions of MATLAB, use line below instead
%rand( 'seed',0);

h = par.h; 
D = abs(min(colrad-abs(X(:,1)-cB(1)),rowrad-abs(X(:,2)-cB(2))));
D2 = D;         
width = lay*h; 
rad = max(D2); 
rel_w = width/rad; 
exp = log2(rel_w); 
sz = -floor(exp)+1; 
prob = 2.^(-lambda*(sz-1:-1:0));

B = D2 < width; 
Ld = min(-ceil(log2(D2./rad)),sz-2)+1;
Ld(B) = sz; 

if sz>0 
P = prob(Ld); 
R = rand(size(P));    
B = R<P; 
end   

Js = 1:size(X,1); 
Xr = X(B,:); 

end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T = LOCAL_Interpolation_Operators(cB,lev,maxlev,params,Xp,X_src,J,k,q)

n = length(X_src); 
if n>k   

form = params.INTERform; 
k_cut = params.k_cut; 
INTpar = params; INTpar.X_source = []; 
INTpar.layers = params.layers+2; 

if params.dim == 3
    Z_src = shape_3D(X_src,params,params.type);      
    X_sk = Z_src(J(1:k),:); 
    X_rs = Z_src(J(k+1:end),:); 
else
    X_sk = X_src(J(1:k),:); 
    X_rs = X_src(J(k+1:end),:);    
end 
              
% INTERPOLATION OPERATORS (VIA LEAST SQUARES) 
% 'dense' - solves T = K[X_ext,X_sk]\K[X_ext,X_rs] densely 
if strcmp(form,'dense') == 1 
    % Dense matrix interpolation operator
    T.U = interpolation_operator(cB,lev,INTpar,Xp,X_sk,X_rs,'box_2_proxy',0); 
    T.V = []; 
else
    % 'lowrank' - (1) [K12.U,K12.V] = LRID(K[X_ext,X_rs]); 
    %             (2) U = K[X_ext,X_sk]\K12.U; V =K12.V;  
    if strcmp(form,'lowrank')                
        if lev < maxlev-1
            % ID-based compression 
            if k < k_cut
                 % Lowrank matrix interpolation operator via LS and ID
                T = interpolation_operator(cB,lev,INTpar,Xp,X_sk,X_rs,'box_2_proxy',1,q); 
            else
                %display(k)  
                INTpar.tau = 1;   %Weighted QR parameter tau
                INTpar.mu = 1e-8; %Tikhonov regularization parameter        
                % FAST (HSS) Lowrank matrix interpolation operator via HSS LS and RandID
                T = interpolation_operator_HSS(cB,lev,INTpar,Xp,X_sk,X_rs,'box_2_proxy',q); 
            end
        else
            % For leaf boxes and their immediate parents, interfaces are
            % too small (or non-existent) 
            T.U = interpolation_operator(cB,lev,INTpar,Xp,X_sk,X_rs,'box_2_proxy',0); 
            T.V = []; 
        end
    end
end
else
   T.U = zeros(k,0); 
   T.V = zeros(0,0); 
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Tf,Jf,kf] = LOCAL_ID_rand_update(T0,J0,k0,X_ext,X,par)

% Parameters
layers = par.layers; 

if par.dim == 3
   X = shape_3D(X,par,par.type);   
end

% Given our initial ID, A_rs ~ A_sk*T. We compute a matvec for the
% residual, Res = A_rs - A_sk*T, of low rank ~ (k-k0), and setup a
% randomized ID: 
J0_sk = J0(1:k0); J0_rs = J0(k0+1:end); 
A_rs = Kernel_Eval(X_ext,X(J0_rs,:),par); 
A_sk = Kernel_Eval(X_ext,X(J0_sk,:),par); 

if strcmp(par.INTERform,'dense') || isempty(T0.V)  
    Res_apply  = @(x) A_rs*x - A_sk*(T0.U*x); 
    Rest_apply = @(x) A_rs.'*x - (T0.U).'*(A_sk.'*x);
else
    Res_apply  = @(x) A_rs*x - A_sk*(T0.U*(T0.V.'*x)); 
    Rest_apply = @(x) A_rs.'*x - T0.V*((T0.U).'*(A_sk.'*x)); 
end

[m,n] = size(A_rs); 
q0 = 8*layers*ceil(log2(k0)); C = q0/2; 

% Rand ID
%fprintf('\n Rand ID extra pts\n')     
[T1,J1] = ID_rand(m,n,Res_apply,Rest_apply,'norm',par.acc,q0,C);

k1 = size(T1,1); 
J1_sk = J1(1:k1); J1_rs = J1(k1+1:end); 

% Update {T,J,k}
kf = k0+k1; 
Jf = [J0_sk ; J0_rs(J1)]; 

if strcmp(par.INTERform,'dense') || isempty(T0.V)    
    Tf.U = [ T0.U(:,J1_rs) - T0.U(:,J1_sk)*T1 ; T1];
    Tf.V = eye(size(Tf.U,2));
else
    % Write Tf as lowrank, may have to recompress. 
    q = size(T0.U,2); 
    Tf.U = [ T0.U zeros(k0,k1); zeros(k1,q) eye(k1) ]; 
    Tf.V = [ (T0.V(J1_rs,:).' - T0.V(J1_sk,:).'*T1) ; T1].';
    
    %Option to Recompress (have to use rand ID instead for optimal complexity)
    %{
    T = LRID(Tf.V.',par.acc);
    Tf.U = Tf.U*T.U; 
    Tf.V = T.V; 
    %}  
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [I_src,count] = Build_I_src(X,leaves,cent_int,N,n2,bisec_rule)
% This function builds the indices of points for each leaf box

nboxes = size(cent_int,1); 
I_src = cell(nboxes,1); 
count = zeros(nboxes,1); 
    
if strcmp(bisec_rule,'uniform')
    % For uniform trees, all leaves are at the same level (depth)
    lev = leaves(1,2); scl = floor(lev/2); numlev = size(leaves,1);
    
    if mod(lev,2) == 0
        % level is even, boxes are square, scl = lev/2;
        coord = ceil((2^(scl-1))*(X + 1));
    else
        % level is odd, boxes are rectangles, scl = (lev-1)/2;
        coord(:,1) = ceil((2^(scl))*(X(:,1) + 1));
        coord(:,2) = ceil((2^(scl-1))*(X(:,2) + 1));
    end
        
    %correcting for points with coordinates == -1. 
    coord = coord + (coord==0); 
    
    box = leaves(:,1).'; 
    tag     = coord(:,1) + n2*coord(:,2) - n2; 
    tag_box = cent_int(box,1) + n2*cent_int(box,2) - n2; 
    ind = repmat(box(tag),numlev,1)==repmat(box(tag_box).',1,N);
    Idx = (1:N).';
    
    count(box) = (sum(ind.')).';
    
    for i=1:numlev
       I_src{box(i)} = Idx(ind(i,:));  
    end
    
else
    % We need to loop over possible levels and build I_src for the
    % corresponding boxes. 
    lev_mn = min(leaves(:,2)); 
    lev_mx = max(leaves(:,2)); 
    Idx = (1:N).';
    
    for lev = lev_mn:lev_mx     
        scl = floor(lev/2);         
        if mod(lev,2) == 0
            % level is even, boxes are square, scl = lev/2;
            coord = ceil((2^(scl-1))*(X + 1));
        else
            % level is odd, boxes are rectangles, scl = (lev-1)/2;
            coord(:,1) = ceil((2^(scl))*(X(:,1) + 1));
            coord(:,2) = ceil((2^(scl-1))*(X(:,2) + 1));
        end
        
        %correcting for points with coordinates == -1. 
        coord = coord + (coord==0); 
    
        box = leaves(leaves(:,2)==lev,1).'; 
        numlev = sum(leaves(:,2)==lev);
        
        tag     = coord(:,1) + n2(lev+1)*coord(:,2) - n2(lev+1); 
        tag_box = cent_int(box,1) + n2(lev+1)*cent_int(box,2) - n2(lev+1); 
        ind = repmat(box(tag),numlev,1)==repmat(box(tag_box).',1,N);
        count(box) = (sum(ind.')).';
    
        for i=1:numlev
            I_src{box(i)} = Idx(ind(i,:));  
        end
    
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I_source,count] = Build_I_source_unif(Xi,N,nbox,nmx,depth,S)
    ind = 1:nmx:N; 
    nX = size(Xi,1);
    lev = depth; 
    scl = floor(lev/2); 
    
    if mod(lev,2) == 0
        % level is even, boxes are square, scl = lev/2;
        coord = ceil((2^(scl-1))*(Xi(:,1:2) + 1));
        %correcting for points with coordinates == -1. 
        coord = coord + (coord==0); 
    else
        % level is odd, boxes are rectangles, scl = (lev-1)/2;
        coord(:,1) = ceil((2^(scl))*(Xi(:,1) + 1));
        coord(:,2) = ceil((2^(scl-1))*(Xi(:,2) + 1));
        %correcting for points with coordinates == -1. 
        coord = coord + (coord==0); 
    end
    
    I_source = cell(nbox,1);
    count = zeros(nbox,1); 
    
    for i=1:nX
        Nbox = S(coord(i,1),coord(i,2)); 
        I_source{Nbox} = ind(i):ind(i)+nmx-1; 
        count(Nbox) = nmx; 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
function [cent,cent_int,C,P,box_numbers,num_bisec,numnodes,numlev,levbox] = bisect_node(centB,numbox,levbox,lev,cent,cent_int,C,P,box_numbers,num_bisec,numnodes,numlev)
    % This function performs the bisection of box B (with center centB, at
    % level lev with coordinates (iB,jB)). 
    
    nchild = 2;
    if mod(levbox(numbox),2) == 0
        % Vertical Bisection
        cta = [-0.5 0 ; 0.5 0];
        scl = lev/2;
    else
        % Horizontal Bisection
        cta = [0 0.5 ; 0 -0.5];
        scl = (lev-1)/2;
    end
    
    
    % Update to data structures (perform the bisection). We assign the last
    % 2 numbers to the new children boxes. 
    % Centers
    cent((1:2) + numnodes*ones(1,2),:) = ones(nchild,1)*centB+cta*(1/2)^(scl);
    % Levels
    levbox((1:2) + numnodes*ones(1,2)) = (lev+1)*ones(1,nchild);
    % Now our box has nchild children
    C(numbox,1:nchild) = (1:2) + numnodes*ones(1,2);
    % These children are, as of now, leaves (no children of their own)
    C((1:2) + numnodes*ones(1,2),:) = zeros(2,2); 
    % Box Numbers
    box_numbers(lev+2,(1:2)+numlev(lev+2)*ones(1,2)) = (1:2) + numnodes*ones(1,2);
    % nchild more on this level
    numlev(lev+2) = numlev(lev+2) + nchild;
    % Parent info for children
    P((1:2) + numnodes*ones(1,2)) = (numbox)*ones(1,nchild);
    
    % We add the children's integer coordinates to S(i,j)
    for ch=1:nchild
        if mod(levbox(numbox)+1,2) == 0
            iC = ceil((2^(scl))*(cent(numnodes+ch,1)+1));
            jC = ceil((2^(scl))*(cent(numnodes+ch,2)+1));
        else
            iC = (ceil((2^(scl+1))*(cent(numnodes+ch,1)+1))+1)/2;
            jC = (ceil((2^(scl))*(cent(numnodes+ch,2)+1))+1)/2;
        end
        
        cent_int(ch + numnodes,:) = [iC jC];     
    end
    
    % Finally, the total number of nodes increases by nchild, and the
    % number of bisections by one. 
    numnodes = numnodes + nchild;
    num_bisec(lev+1) = num_bisec(lev+1) + 1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cent,cent_int,C,P,box_numbers,num_bisec,numnodes,numlev,levbox] = bisect_level(centB,numbox,levbox,lev,cent,cent_int,C,P,box_numbers,num_bisec,numnodes,numlev)
    % This function performs the bisection of box B (with center centB, at
    % level lev with coordinates (iB,jB)). 
    
    nchild = 2^(lev+1);
    scl = floor(lev/2); 
    
    if mod(lev,2) == 0
        % Vertical Bisection
        cta = [-0.5 0 ; 0.5 0];
    else
        % Horizontal Bisection
        cta = [0 0.5 ; 0 -0.5];
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
    ch = 1:nchild; 
    
    if mod(lev+1,2) == 0
        iC = ceil((2^(scl))*(cent(numnodes+ch,1)+1));
        jC = ceil((2^(scl))*(cent(numnodes+ch,2)+1));
    else
        iC = (ceil((2^(scl+1))*(cent(numnodes+ch,1)+1))+1)/2;
        jC = (ceil((2^(scl))*(cent(numnodes+ch,2)+1))+1)/2;
    end
    
    cent_int(ch+numnodes,:) = [iC jC];
    
    % Finally, the total number of nodes increases by nchild, and the
    % number of bisections by one. 
    numnodes = numnodes + nchild;
    num_bisec(lev+1) = 2^lev;

end
