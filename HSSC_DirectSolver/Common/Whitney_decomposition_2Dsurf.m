function [wdec,B] = Whitney_decomposition_2Dsurf(G,cB,rad,layers,h)
%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%

% It is assumed G is a 2D surface, and G \subset [-1,1]^3
N = size(G,1); 
Idx = 1:N; 

% Integer coordinates of curve G on uniform grid 
coord  = [round((-cB(1) + (rad(1)-h/2) + G(:,1))/h) ...
          round((-cB(2) + (rad(2)-h/2) + G(:,2))/h) ...
          round((-cB(3) + (rad(3)-h/2) + G(:,3))/h)]; 

% maximum coordinates
max1 = max(coord(:,1)); 
max2 = max(coord(:,2)); 
max3 = max(coord(:,3)); 
dmax = max([max1 max2 max3]); 

% depth of the decomposition
depth = ceil(log2(dmax));

% Allocation and top level of decomposition
wdec = cell(9,4^depth);
wdec{1,1} = 0; 
wdec{2,1} = NaN; 
wdec{6,1} = cB;
wdec{7,1} = true(N,1); 
wdec{8,1} = N; 
wdec{9,1} = Idx; 

B = false(N,1);  
numnodes = 1; 
numlev = zeros(depth+1,1); 
numlev(1) = 1; 

% Frequency parameter depending on layers
freq = ceil(layers/2); 
sep = 3;  
if size(G,1)>1000 
display(sep)
end 

for lev=0:depth-1
    for box = numnodes-numlev(lev+1)+1:numnodes
        if wdec{8,box}>0
            
            % Box integer coordinates
            box_coord = coord(wdec{9,box},:); 
            [wdec{4,box},id_mn] = min(box_coord); 
            [wdec{5,box},id_mx] = max(box_coord); 
            
            % distance to the boundary
            cent = wdec{6,box};
            dist = min(rad-abs(cent-cB)); 
            boxrad = min(2^(-lev)*rad); 
            
            % If the box contains points and is not well-separated, refine.
            if wdec{8,box}>1 & dist<sep*boxrad
                
                %ind = [id_mn id_mx]; 
                ind = ones(1,0); 
                for c = 1:3
                    not_c = (1:3); not_c(c) = [];
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
                
                % int coordinates for box center
                mid  = [round((-cB(1) + (rad(1)-h/2) + cent(1))/h) ...
                        round((-cB(2) + (rad(2)-h/2) + cent(2))/h) ...
                        round((-cB(3) + (rad(3)-h/2) + cent(3))/h)];
            
                % Split box into 8 
                ch = numnodes+(1:8); 
                cta = [-0.5 -0.5 -0.5; -0.5 0.5 -0.5; 0.5 -0.5 -0.5; 0.5 0.5 -0.5; ...
                       -0.5 -0.5  0.5; -0.5 0.5  0.5; 0.5 -0.5  0.5; 0.5 0.5  0.5]; 
                cta = repmat(rad,8,1).*cta; 
                
                ch_cent = repmat(cent,8,1) + 2^(-lev)*cta;
                
                % boolean vectors for children
                wdec{7,ch(1)} = box_coord(:,1)<mid(1)  & box_coord(:,2)<mid(2)  & box_coord(:,3)<mid(3); 
                wdec{7,ch(2)} = box_coord(:,1)<mid(1)  & box_coord(:,2)>=mid(2) & box_coord(:,3)<mid(3); 
                wdec{7,ch(3)} = box_coord(:,1)>=mid(1) & box_coord(:,2)<mid(2)  & box_coord(:,3)<mid(3);  
                wdec{7,ch(4)} = box_coord(:,1)>=mid(1) & box_coord(:,2)>=mid(2) & box_coord(:,3)<mid(3); 
                wdec{7,ch(5)} = box_coord(:,1)<mid(1)  & box_coord(:,2)<mid(2)  & box_coord(:,3)>=mid(3); 
                wdec{7,ch(6)} = box_coord(:,1)<mid(1)  & box_coord(:,2)>=mid(2) & box_coord(:,3)>=mid(3); 
                wdec{7,ch(7)} = box_coord(:,1)>=mid(1) & box_coord(:,2)<mid(2)  & box_coord(:,3)>=mid(3);  
                wdec{7,ch(8)} = box_coord(:,1)>=mid(1) & box_coord(:,2)>=mid(2) & box_coord(:,3)>=mid(3); 
            
                wdec{3,box} = ch; 
            
                % Children data
                for j=1:8
                    wdec{1,ch(j)} = lev+1; 
                    wdec{2,ch(j)} = box; 
                    wdec{8,ch(j)} = sum(wdec{7,ch(j)}); 
                    wdec{9,ch(j)} = wdec{9,box}(wdec{7,ch(j)}); 
                    wdec{6,ch(j)} = ch_cent(j,:); 
                end
            
                numlev(lev+2) = numlev(lev+2)+8; 
                numnodes = numnodes+8; 
                
                % Update boolean vector of decomposition points
                J = wdec{9,box};   
                B(J(ind)) = true;  
            end
            
            
        end
    end  
end

if numnodes<size(wdec,2)
    wdec = wdec(:,1:numnodes); 
end

maxcrd  = round( 2*(rad-h/2)/h );

% Enforce that all points in boundary layers are included
bdry = (coord(:,1)<layers | maxcrd(1)-coord(:,1)<layers) |...
       (coord(:,2)<layers | maxcrd(2)-coord(:,2)<layers) | ...
       (coord(:,3)<layers | maxcrd(3)-coord(:,3)<layers);
   
    
B = B | bdry;       