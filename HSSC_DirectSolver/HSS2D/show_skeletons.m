function show_skeletons(INFO_TREE,params,lev)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%     FUNCTION CALL:
%         show_skeletons(INFO_TREE,params,lev)
%
%     DESCRIPTION:
%         For a given tree level, this function plots leaf boundaries in red, skeleton points
%         (I_skup) in black, and residual points (I_rsup) in blue (along the interface).
%
%     INPUT:
%         INFO_TREE   <struct>    (output of HSS2D_build_tree_fsym) HSS binary tree.
%         params      <struct>    Kernel parameters (in particular params.X_source is needed).
%         level       <int>       Level in the tree.
%




vlev = INFO_TREE.box_numbers(lev+1,:);
vlev = vlev(vlev>0);

L = length(vlev);
X_source = params.X_source; 

figure; 
hold on
show_tree(INFO_TREE)

for i=1:L
Xsu = X_source(INFO_TREE.BOX(vlev(i)).I_src,:);
Xu = X_source(INFO_TREE.BOX(vlev(i)).I_sk,:);
plot(Xsu(:,1),Xsu(:,2),'ob','MarkerFaceColor','b','MarkerSize',5)
hold on
plot(Xu(:,1),Xu(:,2),'ok','MarkerFaceColor','k','MarkerSize',4)
end
hold off; 

axis equal; axis tight; 

%{
figure

vlevP = INFO_TREE.box_numbers(lev,:);
vlevP = vlevP(vlevP>0);
LP = length(vlevP);
for j=1:LP
     NB = vlevP(j);
     C1 = INFO_TREE.BOX(NB).child(1); 
     C2 = INFO_TREE.BOX(NB).child(2); 
    
     Xu = X_source(INFO_TREE.BOX(NB).I_skup,:); 
     Xu1 = X_source(INFO_TREE.BOX(C1).I_skup,:);
     Xu2 = X_source(INFO_TREE.BOX(C2).I_skup,:);
     
     hold on
     
     plot(Xu1(:,1),Xu1(:,2),'ok','MarkerFaceColor','g','MarkerSize',10)
     plot(Xu2(:,1),Xu2(:,2),'ok','MarkerFaceColor','k','MarkerSize',10)
     plot(Xu(:,1),Xu(:,2),'o','MarkerFaceColor','b','MarkerSize',10)
     
end
%}