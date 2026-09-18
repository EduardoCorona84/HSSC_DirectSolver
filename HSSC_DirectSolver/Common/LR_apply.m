function Y = LR_apply(T,x)
%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%

Y = T.V.'*x; 
Y = T.U*Y; 

return