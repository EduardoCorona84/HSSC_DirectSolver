function [wdec,B] = Whitney_decomposition(X,cB,colrad,rowrad,layers,h,Ng)
%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%

% Integer coordinates on uniform grid 
coord  = [round((-cB(1) + (colrad-h/2) + X(:,1))/h) ...
          round((-cB(2) + (rowrad-h/2) + X(:,2))/h)]; 

% maximum coordinates
max1 = max(coord(:,1)); 
max2 = max(coord(:,2)); 
dmax = max(max1,max2); 

% depth of the decomposition
depth = ceil(log2(dmax))+1;

wdec = cell(depth,1); 

% We determine the first level to be of full boundary layers
wdec{1} = (coord(:,1)<layers | max1-coord(:,1)<layers) |...
          (coord(:,2)<layers | max2-coord(:,2)<layers);
B = wdec{1}; 

% Frequency parameter depending on layers
freq = ceil(layers/2); 

lNg = ceil(log2(Ng)); 

for i=2:depth
    % fixed coordinate ci and modulus mi (either powers of 2 or Ng*(powers
    % of 2) to coincide with grid of Ng x Ng panels. 
    if 2^(i-1)<Ng
        ci = 2^(i-1); mi = Ng; 
    else
        ci = Ng*2^(i-1-lNg); mi = 2*ci;
    end
    
    exp = max(i-freq,0); 
    
    % Add points where each coordinate is fixed 
    wdec{i} = (coord(:,1)==ci | coord(:,1)==max1-ci) ...
        & ( mod(mod(coord(:,2),mi),2^exp)==0 | coord(:,2)==max2-ci);
    
    wdec{i} = wdec{i} | ( (coord(:,2)==ci | coord(:,2)==max2-ci) ...
        & (mod(mod(coord(:,1),mi),2^exp)==0 | coord(:,1)==max1-ci));
    
    B = B | wdec{i}; 
end

display(depth)
display(length(B))
display(sum(B))