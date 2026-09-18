function Z = shape_3D(X,params,type)
%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%


if strcmp(type,'torus')
    r = params.r; 
    w = params.w; 
    m = params.m; 
    b = params.b;
    
    ttc(:,2) =   pi*(X(:,2)+1); 
    ttc(:,1) = m*pi*(X(:,1)+1)+b;        
    Z   = [cos(ttc(:,1)).*(1 + r*cos(ttc(:,2))) ...
            sin(ttc(:,1)).*(1 + r*cos(ttc(:,2))).*(1 + 0.1*cos(ttc(:,1)).*cos(ttc(:,2))) ...
            r*sin(ttc(:,2)).*(1 + 0.1*cos(w*ttc(:,1))) ];       
    
elseif strcmp(type,'patch')
    Z(:,1:2) = X(:,1:2);    
    Z(:,3)   = fnval(params.pp,X(:,1:2).').'; 
end

end