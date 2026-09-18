function C_HSS = HSS1D_updateLR_fsym(A_HSS,U,V,acc)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        C_HSS = HSS1D_updateLR_fsym(A_HSS,U,V,acc)
%
%    DESCRIPTION:
%        This function computes a low rank update UV' to HSS matrix A.
%
%    INPUT:
%        A_HSS is the {50 x nbox} cell array that encodes the HSS structure of matrix A.
%        U,V are the (n x k) factors of the low rank update.
%        acc (double) is the desired accuracy.
%
%    OUTPUT:
%        C_HSS is the {50 x nbox} cell array that encodes the HSS structure of C = A + UV'.
%


global BOX 

nboxes = size(A_HSS,2);
C_HSS = cell(50,nboxes); 

% First step: copy basic tree structure from A_HSS
C_HSS(BOX.DATA,:) = A_HSS(BOX.DATA,:);

% Second pass: Skeleton information and Self/Sibling Interactions
for ibox = nboxes:(-1):1
    if ( (A_HSS{BOX.C1,ibox}<=0) && (A_HSS{BOX.C2,ibox}<=0) ) % ibox has no sons.
        %indices
        indskel = A_HSS{BOX.END1,ibox} - 1 + (1:A_HSS{BOX.END2,ibox});
        
        % Self interaction matrices are just the sum of D^{A} and D^{B}
        C_HSS{BOX.M_SELF,ibox} = A_HSS{BOX.M_SELF,ibox} + U(indskel,:)*(V(indskel,:).'); 
        
        % Concatenate and Re-compress interpolation operators
        % Right / Outgoing Interpolation Matrices
        kA = A_HSS{BOX.KSKEL,ibox}; 
        RA = [eye(kA) A_HSS{BOX.T_UP,ibox}]; RA(:,A_HSS{BOX.J_UP,ibox}) = RA; 
        RB = V(indskel,:).'; 
        RC = [RA ; RB]; 
        
        [T,J] = ID(RC,acc);
        
        % Skeleton size k
        k = size(T,1);
        
        % Record the outgoing skeletons:
        C_HSS{BOX.I_SKUP,ibox} = indskel(J(1:k));
        C_HSS{BOX.T_UP,ibox} = T;
        C_HSS{BOX.J_UP,ibox} = J;
        
        C_HSS{BOX.NSKEL,ibox} = length(indskel);
        C_HSS{BOX.KSKEL,ibox} = k;
        
        % Record pieces of R and L that need to be merged on the next level
        % rho = R^C(:,J(1:k))
        C_HSS{BOX.RHO,ibox} = RC(:,J(1:k));
    else
        ison1       = C_HSS{BOX.C1,ibox};
        ison2       = C_HSS{BOX.C2,ibox};
        if (ison1>0 && ison2>0)
            indskel = [C_HSS{BOX.I_SKUP,ison1},C_HSS{BOX.I_SKUP,ison2}];
            
             % Sibling Interaction Matrices
            kA1 = A_HSS{BOX.KSKEL,ison1}; kA2 = A_HSS{BOX.KSKEL,ison2}; 
            C_HSS{BOX.M_SIB,ison1} = C_HSS{BOX.RHO,ison1}(1:kA1,:).'*A_HSS{BOX.M_SIB,ison1}*C_HSS{BOX.RHO,ison2}(1:kA2,:) +...
                            U(C_HSS{BOX.I_SKUP,ison1},:)*(V(C_HSS{BOX.I_SKUP,ison2},:).');
            C_HSS{BOX.M_SIB,ison2} = C_HSS{BOX.M_SIB,ison1}.'; 
        else
            ison = max([ison1,ison2]); 
            kA1 = A_HSS{BOX.KSKEL,ison}; indskel = C_HSS{BOX.I_SKUP,ison};
        end
                        
        if ibox>1                
        % Concatenate and Re-compress interpolation operators, merging the
        % necessary pieces
        % Right / Outgoing Interpolation Matrices 
        kA = A_HSS{BOX.KSKEL,ibox}; 
        RA = [eye(kA) A_HSS{BOX.T_UP,ibox}]; RA(:,A_HSS{BOX.J_UP,ibox}) = RA; 
        
        if (ison1>0 && ison2>0)
        RC = [RA*[C_HSS{BOX.RHO,ison1}(1:kA1,:) zeros(kA1,C_HSS{BOX.KSKEL,ison2}) ; zeros(kA2,C_HSS{BOX.KSKEL,ison1}) C_HSS{BOX.RHO,ison2}(1:kA2,:)] ;...
              [C_HSS{BOX.RHO,ison1}(kA1+1:end,:) C_HSS{BOX.RHO,ison2}(kA2+1:end,:)]]; 
        else
            RC = [RA*C_HSS{BOX.RHO,ison}(1:kA1,:) ; C_HSS{BOX.RHO,ison}(kA1+1:end,:)]; 
        end
          
        [T,J] = ID(RC,acc);
        
        % Skeleton size k
        k = size(T,1);
        
        % Record the outgoing skeletons:
        C_HSS{BOX.I_SKUP,ibox} = indskel(J(1:k));
        C_HSS{BOX.T_UP,ibox} = T;
        C_HSS{BOX.J_UP,ibox} = J;
        
        C_HSS{BOX.NSKEL,ibox} = length(indskel);
        C_HSS{BOX.KSKEL,ibox} = k; 
        
        % Record pieces of R and L that need to be merged on the next level
        % rho = R^C(:,Jup(1:k))
        C_HSS{BOX.RHO,ibox} = RC(:,J(1:k));
        end
    end
end