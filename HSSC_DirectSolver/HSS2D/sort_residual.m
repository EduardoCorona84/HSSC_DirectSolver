function [X_sort,J] = sort_residual(X,lev)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%     FUNCTION CALL:
%         [X_sort,J] = sort_residual(X,lev)
%
%     DESCRIPTION:
%         This function sorts points along a box interface. If lev is even, the interface is
%         vertical, if it is odd, horizontal. It is generally called on the residual points of a
%         box, which are the points along the shared edge of the children. It is used so that
%         the proxy-residual interaction matrix can be compressed as HSS1D.
%
%     INPUT:
%         X       <(n-k)x2 float>     Residual points that we want sorted.
%         lev     <int>               Level in the tree.
%
%     OUTPUT:
%         X_sort  <(n-k)x2 float>     Array of sorted points, bottom-top or left-right.
%         J       <(n-k)x1 int>       Permutation vector such that X_sort = X(J,:)
%




if mod(lev,2)==0
    [~,J] = sort(-X(:,2));  
else
    [~,J] = sort(-X(:,1));
end

 X_sort = X(J,:);