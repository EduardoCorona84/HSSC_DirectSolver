function [wdec,B] = Whitney_decomposition_1Dcurve(G,C,R,cB,rB,layers,h)
%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%

% It is assumed G is a 1D curve, and G \subset [-1,1]^2
N = size(G,1); 
Idx = 1:N; 

% Integer coordinates of curve G on uniform grid 
coord  = [round((-C(1) + (R(1)-h/2) + G(:,1))/h) ...
          round((-C(2) + (R(2)-h/2) + G(:,2))/h)]; 

% maximum coordinates
max1 = max(coord(:,1)); min1 = min(coord(:,1)); 
max2 = max(coord(:,2)); min2 = min(coord(:,2)); 
dmax = max(max1-min1,max2-min2); 

% depth of the decomposition
depth = ceil(log2(dmax));

% Allocation and top level of decomposition
wdec = cell(9,2^depth);
wdec{1,1} = 0; 
wdec{2,1} = NaN; 
wdec{6,1} = C;
wdec{7,1} = true(N,1); 
wdec{8,1} = N; 
wdec{9,1} = Idx; 

B = false(N,1);   
numnodes = 1; 
numlev = zeros(depth+1,1); 
numlev(1) = 1; 

% Frequency parameter depending on layers
freq = min(ceil(layers/2),2);  

for lev=0:depth-1
    for box = numnodes-numlev(lev+1)+1:numnodes
        if wdec{8,box}>0
            
            box_coord = coord(wdec{9,box},:); 
            [wdec{4,box},id_mn] = min(box_coord); 
            [wdec{5,box},id_mx] = max(box_coord); 
            
            % distance to the boundary  
            cent = wdec{6,box};
            [inf_nrm,i] = max(abs(cent-cB)); 
            dist = abs(rB(i)-inf_nrm); 
            rad = min(2^(-lev)*R); 
            
            % If the box contains points and is not well-separated, refine.
            if wdec{8,box}>1 & dist<freq*1.5*rad
                ind = [id_mn id_mx];   
                
                % int coordinates for box center
                cBc  = [round((-C(1) + (R(1)-h/2) + cent(1))/h) ...
                        round((-C(2) + (R(2)-h/2) + cent(2))/h)];
                rowmid = cBc(1); colmid = cBc(2); 
                %rowmid = round((wdec{5,box}(1)+wdec{4,box}(1))/2);
                %colmid = round((wdec{5,box}(2)+wdec{4,box}(2))/2);
            
                % Split box into 4 if box is well-separated
                ch = numnodes+(1:4); 
                cta = [-0.5 -0.5; -0.5 0.5; 0.5 -0.5; 0.5 0.5]; 
                ch_cent = repmat(cent,4,1) + 2^(-lev)*cta;
                
                wdec{7,ch(1)} = (box_coord(:,1)<rowmid & box_coord(:,2)<colmid); 
                wdec{7,ch(2)} = (box_coord(:,1)<rowmid & box_coord(:,2)>=colmid); 
                wdec{7,ch(3)} = (box_coord(:,1)>=rowmid & box_coord(:,2)<colmid);  
                wdec{7,ch(4)} = (box_coord(:,1)>=rowmid & box_coord(:,2)>=colmid); 
            
                wdec{3,box} = ch; 
            
                for j=1:4
                    wdec{1,ch(j)} = lev+1; 
                    wdec{2,ch(j)} = box; 
                    wdec{8,ch(j)} = sum(wdec{7,ch(j)}); 
                    wdec{9,ch(j)} = wdec{9,box}(wdec{7,ch(j)}); 
                    wdec{6,ch(j)} = ch_cent(j,:); 
                end
            
                numlev(lev+2) = numlev(lev+2)+4; 
                numnodes = numnodes+4; 
                
                %{
                ind = ones(1,0); 
                for c = 1:2
                    not_c = 3-c; 
                    sn = 1:length(wdec{9,box}); 
                    
                    smn = sn(box_coord(:,c)==wdec{4,box}(c)); 
                    
                    if length(smn)>1
                        [~,id_mn] = min(box_coord(box_coord(:,c)==wdec{4,box}(c),not_c)); 
                        id_mn = smn(id_mn); 
                    else
                        id_mn = smn; 
                    end
                    
                    smx = sn(box_coord(:,c)==wdec{5,box}(c)); 
                    
                    if length(smx)>1
                        [~,id_mx] = max(box_coord(box_coord(:,c)==wdec{5,box}(c),not_c)); 
                        id_mx = smx(id_mx); 
                    else
                        id_mx = smx; 
                    end
                    
                    ind = [ind id_mn id_mx]; 
                end
                %}
                
                J = wdec{9,box}; 
                B(J(ind)) = true;   
            end
            
            
        end
    end  
end

% Enforce that all points in boundary layers are included
inf_nrm = abs(G-repmat(cB,N,1)); 
dist = min(abs(repmat(rB,N,1)-inf_nrm),[],2); 
bdry = dist<layers*h & (inf_nrm(:,1)<rB(1)+layers*h & inf_nrm(:,2)<rB(2)+layers*h);

%bdry = (coord(:,1)-min1<layers | max1-coord(:,1)<layers) |...
%          (coord(:,2)-min1<layers | max2-coord(:,2)<layers);

B = B | bdry;       