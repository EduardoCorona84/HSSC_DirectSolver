function params = HSS_tree_parameters(pot,kh,sym,Ng,np,acc,lay,n_cut,TI)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%     FUNCTION CALL:
%         params = HSS_tree_parameters(pot,kh,sym,Ng,np,acc,lay,n_cut,TI)
%
%     DESCRIPTION:
%         This function lets you input parameters for the kernel, binary tree, and compression.
%         The default discretization is taken to be a regular grid of [-1,1]x[-1,1] consisting
%         np x np panels, each with Ng x Ng tensor gaussian nodes. This can be changed by
%         redefining params.X_source after the function terminates. The return value params is
%         used as an argument in most HSS2D functions, however it must be updated during the
%         BUILD TREE stage (after HSS2D_bintree_* is called).
%
%     INPUT:
%         pot     <string>    (.flagpot) Name of the kernel for matrix entries K[p,q].
%                                 See Kernel_Eval for current options.
%         kh      <int>       (.kh) Wavenumber, if applicable. E.g. when evaluating low
%                                 frequency oscillatory kernels such as the 2D Helmholtz free-
%                                 space Green's function.
%         sym     <bool>      (.sym) Whether or not the kernel is symmetric.
%         Ng      <int>       (.Ng) Number of points along one side of a panel (typically 5 - 9).
%         np      <int>       Number of panels along one side of the domain (for uniform tree).
%         acc     <float>     (.acc) Accuracy of IDs for leaf boxes.
%         lay     <int>       (.layers) Number of boundary layers used for skeleton sets.
%                                 If kh<30, we suggest lay = 1 or 2 if acc is greater
%                                 or less than 1e-7. If kh>30 we suggest lay=3.
%         n_cut   <int>       (.n_cut) Minimum skeleton set size to switch to HSS form.
%         TI      <bool>      (.transinv) Whether or not kernel is translation invariant.
%
%
%     OUTPUT:
%         params  <struct>    Kernel, tree, and compression parameters for HSS2D algorithms.
%                                 Default values are shown below.
%
%         Kernel Parameters
%             .flagpot = pot
%             .kh = kh
%             .sym = sym
%             .transinv = TI
%
%         Discretization Parameters
%             .Ng
%             .h = 2/(Ng * np)    Grid spacing of points.
%             .X_source           Source points in discretization. Calculated as guassian
%                                     quadrature nodes in each panel in [-1,1]x[-1,1].
%             .order = 4;         Order of quadrature.
%             .dr_weights         Duan-Rohklin quadrature weights.
%
%
%         Tree Parameters
%             .max_particles = Ng*Ng      Maximum number of source points allowed in a leaf
%                                             box.
%             .depth ~ 2log(np);          Maximum depth of quadtree so that a leaf box
%                                             equals a panel.
%             .bisec_rule = 'uniform'     {'uniform', 'maxpart'} Binary tree stopping criterion.
%                                         - 'uniform' creates a unifom quadtree of .depth levels
%                                         - 'maxpart' is an adaptive build where a box is
%                                             bisected only if it contains greater than
%                                             max_particles source points and is at a level
%                                             greater than .depth.
%
%         Compression Parameters
%             .acc = acc
%             .layers = lay
%             .n_cut = n_cut
%             .DOskeletons = true     After the tree is constructed with box centers,
%                                         this determines whether to proceed in the
%                                         compression of A by computing skeleton/index sets and
%                                         interpolation operators.
%             .skel_rule = 'bdry'     {'ID', 'wdec', 'brdy'} Method to construct skeleton sets.
%                                     - 'ID' Skeleton found by interpolatory decomposition.
%                                     - 'wdec' Skeleton pre-determined by Whitney decomposition.
%                                     - 'brdy' Skeleton pre-determined to be .layers layers of
%                                         points along the boundary of the box.
%             .skel_extra_pts = 0     Boolean value that specifies whether to add extra points
%                                         to the preset skeleton sets. The extra points are
%                                         computed by a randomized ID and they are used to keep
%                                         interpolation accurate. This parameter should be set
%                                         to true for oscillatory kernels or kernels that don't
%                                         fulfill a Green's identity in the domain.
%             .INTERform = 'lowrank'  {'dense', 'lowrank'} Determines computation and storage of
%                                         interpolation operators.
%             .k_cut = 2500           Minimum skeleton set size to compute L,R matrices with
%                                         HSS acceleration.
%             .proxy_rule = 'rand'    {'wdec', 'rand', 'brdy', neigh'} Method to construct
%                                         proxies for equivalent density acceleration.
%                                     - 'wdec' Whitney decomposition.
%                                     - 'rand' Exterior boundary layers plus random points
%                                         selected from an exponentially decaying distribution.
%                                     - 'brdy' Exterior boundary layers of the box.
%                                     - 'neigh' Proxy taken from skeleton points of neighbors.
%             .lambda = 1             For .proxy_rule='rand', the distribution of proxy points
%                                         ~ 2^-(lambda*r) where r is the distance.
%             .dim = 2                Dimension of the kernel.
%             .proxy = 0              Proxy to skeleton interactions for NTI kernels. Specifies
%                                         perturbation functions b and c.
%             .pp = []                Spline data for 3D examples.
%




hp = 2/np; 
hh = 2/(Ng * np);
ntot    = np*np*Ng*Ng;  % Total number of nodes.
display(ntot)
display(lay)
nlevels = round(2*(log(np)/log(2) - 1)+3);

% Parameters 
params.Ng = Ng; 
params.bisec_rule = 'uniform';
params.depth = nlevels - 1;
params.DOskeletons = 1; 
params.max_particles = Ng*Ng; 
params.skel_rule = 'bdry'; 
params.skel_extra_pts = 0; 
params.proxy_rule = 'bdry';   
params.lambda = 2;  
params.layers = lay; 
params.h = hh;
params.flag_pot = pot;
params.kh = kh; 
params.acc = acc; 
params.n_cut = n_cut; 
params.INTERform = 'lowrank'; 
params.transinv = TI;
params.sym = sym; 
params.k_cut = 2500; 
params.dim = 2;
params.m = 1; 
params.type = 'plane';    
params.proxy = 0; 
params.pp = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build X_source
%%% Construct the Gaussian nodes for a single panel.
z = linspace(0.5*hh,Ng*hh-0.5*hh,Ng);
w = hh*ones(Ng,1);

%%% Construct the tensor product quadrature for a panel.
[ZZ1,ZZ2] = meshgrid(z);
%WW        = w*w';
J         = Ng*Ng;
zzloc     = [reshape(ZZ1,1,J);...
             reshape(ZZ2,1,J)];
%wwloc     = reshape(WW,J,1);

[X1,Y1] = meshgrid(-1:hp:1-hp,-1:hp:1-hp);
Xmin = [X1(:) Y1(:)];

for i=1:np*np
X_source((J*(i-1)+1):J*i,:) = [Xmin(i,1)*ones(J,1) Xmin(i,2)*ones(J,1)] + zzloc';
end

params.X_source = X_source; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Helmholtz kernel, generate Duan-Rokhlin correction weights

% Duan-Rokhlin weights information
if kh > 0 
   params.order = 4; 	% quadrature order, set = 0 to forgo DR weights
   [dr_J,dr_weights,iistencil] = DR_stencil(hh,kh,np*Ng,params.order);
   params.dr_weights = dr_weights; 
   if params.order > 4
      params.dr_J = dr_J; 
      params.iistencil = iistencil; 
   end
end
