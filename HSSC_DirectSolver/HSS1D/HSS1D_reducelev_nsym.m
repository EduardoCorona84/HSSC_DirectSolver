function ENEW = HSS1D_reducelev_nsym(E,acc)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        ENEW = HSS1D_reducelev_nsym(E,acc)
%
%    DESCRIPTION:
%        This function reduces one level of the non-symmetric HSS matrix E by eliminating leaf
%        nodes, updating their parent's interpolation matrix R, and recompressing up the tree.
%
%    INPUT:
%        E is a {50 x nbox} cell array containing HSS info for matrix E.
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
if ( (E{BOX.C1,ibox}<=0) && (E{BOX.C2,ibox}<=0) )
    if ibox>1
    k = E{BOX.KSKEL,ibox}; 
    Jout = E{BOX.J_UP,ibox}; Jin = E{BOX.J_DN,ibox}; 
    R = [eye(k) E{BOX.T_UP,ibox}]; L = [eye(k) E{BOX.T_DN,ibox}].';         
    R(:,Jout) = R; 
    L(Jin,:) = L; 
    
    % Eliminate Leaf Node
    par = E{BOX.PARENT,ibox}; child = [E{BOX.C1,par} E{BOX.C2,par}];  
    sib = child(child~=ibox);
    
    if sib<=0 
        ENEW{BOX.NSKEL,ibox} = 0; 
        ENEW{BOX.RHO,ibox} = R; 
        ENEW{BOX.LAMBDA,ibox} = L; 
    elseif ( (E{BOX.C1,sib}<=0) && (E{BOX.C2,sib}<=0) )
        ENEW{BOX.NSKEL,ibox} = 0; 
        ENEW{BOX.RHO,ibox} = R; 
        ENEW{BOX.LAMBDA,ibox} = L; 
    else
        ENEW(8:40,ibox) = E(8:40,ibox); 
        ENEW{BOX.RHO,ibox} = eye(size(R,1)); 
        ENEW{BOX.LAMBDA,ibox} = eye(size(L,2)); 
    end
    end
else
    if ibox>1
        k = E{BOX.KSKEL,ibox}; 
        Jout = E{BOX.J_UP,ibox}; Jin = E{BOX.J_DN,ibox}; 
        R = [eye(k) E{BOX.T_UP,ibox}]; L = [eye(k) E{BOX.T_DN,ibox}].';         
        R(:,Jout) = R; 
        L(Jin,:) = L;    
        %Non-leaf boxes
        %Two children case 
        if ( (E{BOX.C1,ibox}>0) && (E{BOX.C2,ibox}>0) )
            ison1 = E{BOX.C1,ibox}; ison2 = E{BOX.C2,ibox};
            if (E{BOX.C1,ison1}<=0 && E{BOX.C2,ison1}<=0) && (E{BOX.C1,ison2}<=0 && E{BOX.C2,ison2}<=0)
                indskel = E{BOX.END1,ibox} - 1 + (1:E{BOX.END2,ibox});
                indskel_out = indskel; indskel_in = indskel; 
            elseif (E{BOX.C1,ison1}>0 || E{BOX.C2,ison1}>0)
                ind2 = E{BOX.END1,ison2} - 1 + (1:E{BOX.END2,ison2});
                indskel_out = [ENEW{BOX.I_SKUP,ison1} ind2];
                indskel_in  = [ENEW{BOX.I_SKDN,ison1} ind2];
            elseif (E{BOX.C1,ison2}>0 || E{BOX.C2,ison2}>0)
                ind1 = E{BOX.END1,ison1} - 1 + (1:E{BOX.END2,ison1});
                indskel_out = [ind1 ENEW{BOX.I_SKUP,ison2}];
                indskel_in  = [ind1 ENEW{BOX.I_SKDN,ison2}];
            else
                indskel_out = [ENEW{BOX.I_SKUP,ison1},ENEW{BOX.I_SKUP,ison2}];
                indskel_in  = [ENEW{BOX.I_SKDN,ison1},ENEW{BOX.I_SKDN,ison2}];
            end
        
            R_c1 = ENEW{BOX.RHO,ison1}; k1 = size(R_c1,1); 
            R_c2 = ENEW{BOX.RHO,ison2};
            
            L_c1 = ENEW{BOX.LAMBDA,ison1}; 
            L_c2 = ENEW{BOX.LAMBDA,ison2}; 
        
            R = [R(:,1:k1)*R_c1 R(:,k1+1:end)*R_c2];
            L = [L_c1*L(1:k1,:) ; L_c2*L(k1+1:end,:)];
        else
        %One Child case
            ison = max([E{BOX.C1,ibox} E{BOX.C2,ibox}]);
            if (E{BOX.C1,ison}<=0 && E{BOX.C2,ison}<=0)
                indskel = E{BOX.END1,ibox} - 1 + (1:E{BOX.END2,ibox});
                indskel_out = indskel; indskel_in = indskel; 
            else
                indskel_out = ENEW{BOX.I_SKUP,ison}; 
                indskel_in = ENEW{BOX.I_SKDN,ison}; 
            end
       
            R = R*ENEW{BOX.RHO,ison}; 
            L = ENEW{BOX.LAMBDA,ison}*L; 
        end
    
        %Re-compression
        %acc = 1e-12; 
        
        [Tout,Jout] = ID(R,acc);
        [Tin,Jin] = ID(L.',acc);
        
        k = max(size(Tout,1),size(Tin,1)); 
        
        if ~(k == size(Tout,1))
            [Tout,Jout] = ID(R, k);
        elseif ~(k == size(Tin,1))
            [Tin, Jin ] = ID(L.',k);
        end 
        
        ENEW{BOX.NSKEL,ibox} = length(indskel_out); 
        ENEW{BOX.KSKEL,ibox} = k; 
        
        ENEW{BOX.I_SKUP,ibox} = indskel_out(Jout(1:k)); %I_skup
        ENEW{BOX.T_UP,ibox} = Tout; % Tout(I_skup,I_rsup)
        ENEW{BOX.J_UP,ibox} = Jout; %Jout
        
        ENEW{BOX.I_SKDN,ibox} = indskel_in(Jin(1:k)); %I_skup
        ENEW{BOX.T_DN,ibox} = Tin; % Tout(I_skup,I_rsup)
        ENEW{BOX.J_DN,ibox} = Jin; %Jout
    
        ENEW{BOX.RHO,ibox} = R(:,Jout(1:k)); 
        ENEW{BOX.LAMBDA,ibox} = L(Jin(1:k),:); 
    end
    
    % SELF and SIBLING INFO
    %Two children case 
    if ( (E{BOX.C1,ibox}>0) && (E{BOX.C2,ibox}>0) )
        ison1 = E{BOX.C1,ibox}; ison2 = E{BOX.C2,ibox};
        ENEW{BOX.M_SIB,ison1} = ENEW{BOX.LAMBDA,ison1}*E{BOX.M_SIB,ison1}*ENEW{BOX.RHO,ison2}; 
        ENEW{BOX.M_SIB,ison2} = ENEW{BOX.LAMBDA,ison2}*E{BOX.M_SIB,ison2}*ENEW{BOX.RHO,ison1}; 
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
    end
    if ENEW{BOX.C2,j}>0
        ENEW{BOX.C2,j} = sumfullE(ENEW{BOX.C2,j});
    end
end