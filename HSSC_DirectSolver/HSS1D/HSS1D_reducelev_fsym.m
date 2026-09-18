function ENEW = HSS1D_reducelev_fsym(E,acc)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        ENEW = HSS1D_reducelev_fsym(E,acc)
%
%    DESCRIPTION:
%        This function reduces one level of the full-symmetric HSS matrix E by eliminating leaf
%        nodes, updating their parent's interpolation matrix R, and recompressing up the tree.
%
%    INPUT:
%        E is a {50 x nbox} cell array containing HSS info for matrix E
%        acc (double) is the desired accuracy.
%
%    OUTPUT:
%        E_NEW is the {50 x nbox_NEW} cell array containing HSS info for reduced E.
%


global BOX 

ENEW(BOX.DATA,:) = E(BOX.DATA,:); 
nboxes = size(E,2); 

for ibox=nboxes:-1:1
% Leaf boxes
if ( (E{BOX.C1,ibox}<=0) & (E{BOX.C2,ibox}<=0) )
    if ibox>1
    k = E{BOX.KSKEL,ibox}; 
    J_up = E{BOX.J_UP,ibox};  
    R = [eye(k) E{BOX.T_UP,ibox}];        
    R(:,J_up) = R;    
    % Eliminate Leaf Node
    par = E{BOX.PARENT,ibox}; child = [E{BOX.C1,par} E{BOX.C2,par}];  
    sib = child(child~=ibox);      
    
    if sib<=0 
        ENEW{BOX.NSKEL,ibox} = 0; 
        ENEW{BOX.RHO,ibox} = R; 
    elseif ( (E{BOX.C1,sib}<=0) && (E{BOX.C2,sib}<=0) )
        ENEW{BOX.NSKEL,ibox} = 0; 
        ENEW{BOX.RHO,ibox} = R; 
    else
        ENEW(8:40,ibox) = E(8:40,ibox); 
        ENEW{BOX.RHO,ibox} = eye(size(R,1)); 
    end
    end
else
    if ibox>1
        k = E{BOX.KSKEL,ibox}; 
        J_up = E{BOX.J_UP,ibox}; 
        R = [eye(k) E{BOX.T_UP,ibox}];        
        R(:,J_up) = R;    
        %Non-leaf boxes
        %Two children case 
        if ( (E{BOX.C1,ibox}>0) && (E{BOX.C2,ibox}>0) )
            ison1 = E{BOX.C1,ibox}; ison2 = E{BOX.C2,ibox};
            if (E{BOX.C1,ison1}<=0 && E{BOX.C2,ison1}<=0) && (E{BOX.C1,ison2}<=0 && E{BOX.C2,ison2}<=0)
                indskel = E{BOX.END1,ibox} - 1 + (1:E{BOX.END2,ibox});
            elseif (E{BOX.C1,ison1}>0 || E{BOX.C2,ison1}>0)
                ind2 = E{BOX.END1,ison2} - 1 + (1:E{BOX.END2,ison2});
                indskel = [ENEW{BOX.I_SKUP,ison1} ind2];
            elseif (E{BOX.C1,ison2}>0 || E{BOX.C2,ison2}>0)
                ind1 = E{BOX.END1,ison1} - 1 + (1:E{BOX.END2,ison1});
                indskel = [ind1 ENEW{BOX.I_SKUP,ison2}];
            else
                indskel = [ENEW{BOX.I_SKUP,ison1},ENEW{BOX.I_SKUP,ison2}];
            end
        
            R_c1 = ENEW{BOX.RHO,ison1}; k1 = size(R_c1,1); 
            R_c2 = ENEW{BOX.RHO,ison2};
        
            R = [R(:,1:k1)*R_c1 R(:,k1+1:end)*R_c2];
        else
        %One Child case
            ison = max([E{BOX.C1,ibox} E{BOX.C2,ibox}]);
            if (E{BOX.C1,ison}<=0 && E{BOX.C2,ison}<=0)
                indskel = E{BOX.END1,ibox} - 1 + (1:E{BOX.END2,ibox});
            else
                indskel = ENEW{BOX.I_SKUP,ison};
            end
       
            R = R*ENEW{BOX.RHO,ison}; 
        end
    
        %Re-compression
        %acc = 1e-12; 
        
        [T_up,J_up] = ID(R,acc);
        k = size(T_up,1); 
        
        ENEW{BOX.NSKEL,ibox} = length(indskel); 
        ENEW{BOX.KSKEL,ibox} = k; 
        ENEW{BOX.J_UP,ibox} = J_up; %J_up
        ENEW{BOX.I_SKUP,ibox} = indskel(J_up(1:k)); %I_skup
        ENEW{BOX.T_UP,ibox} = T_up; % T_up(I_skup,I_rsup)
    
        ENEW{BOX.RHO,ibox} = R(:,J_up(1:k)); 
    end
    
    % SELF and SIBLING INFO
    %Two children case 
    if ( (E{BOX.C1,ibox}>0) && (E{BOX.C2,ibox}>0) )
        ison1 = E{BOX.C1,ibox}; ison2 = E{BOX.C2,ibox};  
        ENEW{BOX.M_SIB,ison1} = ENEW{BOX.RHO,ison1}.'*E{BOX.M_SIB,ison1}*ENEW{BOX.RHO,ison2}; 
        ENEW{BOX.M_SIB,ison2} = ENEW{BOX.M_SIB,ison1}.'; 
        if (E{BOX.C1,ison1}<=0 && E{BOX.C2,ison1}<=0) && (E{BOX.C1,ison2}<=0 && E{BOX.C2,ison2}<=0)
        % Case I: Children are leaves
            ENEW{BOX.M_SELF,ibox} = [E{BOX.M_SELF,ison1} ENEW{BOX.M_SIB,ison1} ; ENEW{BOX.M_SIB,ison2} E{BOX.M_SELF,ison2}]; 
            % Now ibox is a leaf
            ENEW{BOX.C1,ibox} = -1; ENEW{BOX.C2,ibox} = -1;     
        end
    else
        ison = max([E{BOX.C1,ibox} E{BOX.C2,ibox}]); 
        if (E{BOX.C1,ison}<=0 && E{BOX.C2,ison}<=0)
            % Case I: Child is a leaf
            ENEW{BOX.M_SELF,ibox} = E{BOX.M_SELF,ison};
            ENEW{BOX.C1,ibox} = -1; ENEW{BOX.C2,ibox} = -1;     
        end
    end
end

end

% Clean up
fullE = zeros(1,nboxes); sumfullE = fullE; 
fullE(1) = 1; sumfullE(1) = 1; 
for j=2:nboxes
    if ENEW{BOX.NSKEL,j} > 0
        fullE(j) = 1;
    end
    sumfullE(j) = sum(fullE(1:j)); 
end

% Eliminate empty boxes
ENEW = ENEW(:,fullE==1); 

%New parent / children info: 
for j=1:size(ENEW,2)

    if ENEW{BOX.C1,j}>0
        ENEW{BOX.C1,j} = sumfullE(ENEW{BOX.C1,j});
        if ENEW{BOX.C1,j}>0
            ENEW{BOX.PARENT,ENEW{BOX.C1,j}} = j; 
        else
            ENEW{BOX.C1,j} = -1; 
        end
    end
    if ENEW{BOX.C2,j}>0
        ENEW{BOX.C2,j} = sumfullE(ENEW{BOX.C2,j});
        if ENEW{BOX.C2,j}>0
            ENEW{BOX.PARENT,ENEW{BOX.C2,j}} = j; 
        else
            ENEW{BOX.C2,j} = -1; 
        end
    end    
end