% Test HSS2D_bintree for the case of a uniform grid on the unit box. We
% check both tree information and skeleton information. 

fprintf('\n Binary 2D Tree Build \n')
tic; 
if params.sym == 0
   [INFO_TREE]=HSS2D_build_tree_nsym(params); 
else
   [INFO_TREE]=HSS2D_build_tree_fsym(params); 
end
T_tree = toc; 

fprintf('\n Tree Build Time O(N): %e \n',T_tree);

%{
% Optional: TREE visualization
for lev = INFO_TREE.depth:-1:0
figure
show_skeletons(INFO_TREE,params,lev);
end
 %}

% Check interpolation matrices
Err = Test_HSS2D_skeletons(INFO_TREE,params);

fprintf('\n Interpolation matrix log10(Errors) by level 0-depth \n')
display(log10(Err))