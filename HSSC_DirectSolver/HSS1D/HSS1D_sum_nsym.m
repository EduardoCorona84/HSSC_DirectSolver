function C_HSS = HSS1D_sum_nsym(A_HSS,B_HSS,acc)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        C_HSS = HSS1D_sum_nsym(A_HSS,B_HSS,acc)
%
%    DESCRIPTION:
%        This function adds two non-symmetric HSS matrices A and B that share source and target
%        trees.
%
%    INPUT:
%        A_HSS is the {50 x nbox} cell array that encodes the HSS structure of matrix A.
%        B_HSS is the {50 x nbox} cell array that encodes the HSS structure of matrix B.
%        acc (double) is the desired accuracy.
%
%    OUTPUT:
%        C_HSS is the {50 x nbox} cell array that encodes the HSS structure of matrix A + B.
%


global BOX 

nboxes = max(size(A_HSS,2),size(B_HSS,2));
C_HSS = cell(50,nboxes); 

% First step: copy basic tree structure from A_HSS / B_HSS
C_HSS(BOX.DATA,:) = A_HSS(BOX.DATA,:);

%TODO: for the case in which A_HSS is deeper than B_HSS (or viceversa),
%we can extend the other trivially bisecting nodes and adding empty
%"skeletons" and interpolation matrices. 


% Second pass: Skeleton information and Self/Sibling Interactions
for ibox = nboxes:(-1):1
    if ( (A_HSS{BOX.C1,ibox}<=0) && (A_HSS{BOX.C2,ibox}<=0) ) % ibox has no sons.
        %indices
        indskel = A_HSS{BOX.END1,ibox} - 1 + (1:A_HSS{BOX.END2,ibox});
        
        % Self interaction matrices are just the sum of D^{A} and D^{B}
        C_HSS{BOX.M_SELF,ibox} = A_HSS{BOX.M_SELF,ibox} + B_HSS{BOX.M_SELF,ibox}; 
        
        % Concatenate and Re-compress interpolation operators
        kA = A_HSS{BOX.KSKEL,ibox}; kB = B_HSS{BOX.KSKEL,ibox}; 
        % Right / Outgoing Interpolation Matrices
        RA = [eye(kA) A_HSS{BOX.T_UP,ibox}]; RA(:,A_HSS{BOX.J_UP,ibox}) = RA; 
        RB = [eye(kB) B_HSS{BOX.T_UP,ibox}]; RB(:,B_HSS{BOX.J_UP,ibox}) = RB; 
        % Left / Outgoing Interpolation Matrices
        LA = [eye(kA) A_HSS{BOX.T_DN,ibox}].'; LA(A_HSS{BOX.J_DN,ibox},:) = LA; 
        LB = [eye(kB) B_HSS{BOX.T_DN,ibox}].'; LB(B_HSS{BOX.J_DN,ibox},:) = LB; 
        
        RC = [RA ; RB]; 
        LC = [LA   LB];
        
        [Tout,Jout] = ID(RC,acc);
        [Tin,Jin] = ID(LC.',acc); 
        
        % Skeleton size k
        k = max(size(Tout,1),size(Tin,1));
        
        if ~(k == size(Tout,1))
            [Tout,Jout] = ID(RC, k);
        elseif ~(k == size(Tin,1))
            [Tin, Jin ] = ID(LC.',k);
        end
        
        % Record the outgoing skeletons:
        C_HSS{BOX.I_SKUP,ibox} = indskel(Jout(1:k));
        C_HSS{BOX.T_UP,ibox} = Tout;
        C_HSS{BOX.J_UP,ibox} = Jout;
        
        % Record the ingoing skeletons
        C_HSS{BOX.I_SKDN,ibox} = indskel(Jin(1:k));
        C_HSS{BOX.T_DN,ibox} = Tin;
        C_HSS{BOX.J_DN,ibox} = Jin;
        
        C_HSS{BOX.NSKEL,ibox} = length(indskel);
        C_HSS{BOX.KSKEL,ibox} = k;
        
        % Record pieces of R and L that need to be merged on the next level
        % rho = R^C(:,J(1:k))
        C_HSS{BOX.RHO,ibox} = RC(:,Jout(1:k));
        C_HSS{BOX.LAMBDA,ibox} = LC(Jin(1:k),:); 
    else
        ison1       = C_HSS{BOX.C1,ibox};
        ison2       = C_HSS{BOX.C2,ibox};
        if (ison1>0 && ison2>0)
            indskel_out = [C_HSS{BOX.I_SKUP,ison1},C_HSS{BOX.I_SKUP,ison2}];
            indskel_in  = [C_HSS{BOX.I_SKDN,ison1},C_HSS{BOX.I_SKDN,ison2}];
            
             % Sibling Interaction Matrices
            kA1 = A_HSS{BOX.KSKEL,ison1}; kB1 = B_HSS{BOX.KSKEL,ison1};
            kA2 = A_HSS{BOX.KSKEL,ison2}; kB2 = B_HSS{BOX.KSKEL,ison2};
            
            C_HSS{BOX.M_SIB,ison1} = C_HSS{BOX.LAMBDA,ison1}(:,1:kA1)*A_HSS{BOX.M_SIB,ison1}*C_HSS{BOX.RHO,ison2}(1:kA2,:) +...
                            C_HSS{BOX.LAMBDA,ison1}(:,kA1+1:end)*B_HSS{BOX.M_SIB,ison1}*C_HSS{BOX.RHO,ison2}(kA2+1:end,:);
            C_HSS{BOX.M_SIB,ison2} = C_HSS{BOX.LAMBDA,ison2}(:,1:kA2)*A_HSS{BOX.M_SIB,ison2}*C_HSS{BOX.RHO,ison1}(1:kA1,:) +...
                            C_HSS{BOX.LAMBDA,ison2}(:,kA2+1:end)*B_HSS{BOX.M_SIB,ison2}*C_HSS{BOX.RHO,ison1}(kA1+1:end,:); 
        else
            ison = max([ison1,ison2]); 
            kA1 = A_HSS{BOX.KSKEL,ison}; kB1 = B_HSS{BOX.KSKEL,ison}; 
            indskel_out = C_HSS{BOX.I_SKUP,ison};
            indskel_in  = C_HSS{BOX.I_SKDN,ison};
        end
                        
        if ibox>1                
        % Concatenate and Re-compress interpolation operators, merging the
        % necessary pieces
        kA = A_HSS{BOX.KSKEL,ibox}; kB = B_HSS{BOX.KSKEL,ibox}; 
        % Right / Outgoing Interpolation Matrices
        RA = [eye(kA) A_HSS{BOX.T_UP,ibox}]; RA(:,A_HSS{BOX.J_UP,ibox}) = RA; 
        RB = [eye(kB) B_HSS{BOX.T_UP,ibox}]; RB(:,B_HSS{BOX.J_UP,ibox}) = RB; 
        % Left / Outgoing Interpolation Matrices
        LA = [eye(kA) A_HSS{BOX.T_DN,ibox}].'; LA(A_HSS{BOX.J_DN,ibox},:) = LA; 
        LB = [eye(kB) B_HSS{BOX.T_DN,ibox}].'; LB(B_HSS{BOX.J_DN,ibox},:) = LB; 
        
        if (ison1>0 && ison2>0)
            RC = [RA*[C_HSS{BOX.RHO,ison1}(1:kA1,:) zeros(kA1,C_HSS{BOX.KSKEL,ison2}) ; zeros(kA2,C_HSS{BOX.KSKEL,ison1}) C_HSS{BOX.RHO,ison2}(1:kA2,:)] ;...
              RB*[C_HSS{BOX.RHO,ison1}(kA1+1:end,:) zeros(kB1,C_HSS{BOX.KSKEL,ison2}) ; zeros(kB2,C_HSS{BOX.KSKEL,ison1}) C_HSS{BOX.RHO,ison2}(kA2+1:end,:)]]; 
          
            LC = [[C_HSS{BOX.LAMBDA,ison1}(:,1:kA1) zeros(C_HSS{BOX.KSKEL,ison1},kA2) ; zeros(C_HSS{BOX.KSKEL,ison2},kA1) C_HSS{BOX.LAMBDA,ison2}(:,1:kA2)]*LA ...
              [C_HSS{BOX.LAMBDA,ison1}(:,kA1+1:end) zeros(C_HSS{BOX.KSKEL,ison1},kB2) ; zeros(C_HSS{BOX.KSKEL,ison2},kB1) C_HSS{BOX.LAMBDA,ison2}(:,kA2+1:end)]*LB]; 
        else
            RC = [RA*C_HSS{BOX.RHO,ison}(1:kA1,:) ; RB*C_HSS{BOX.RHO,ison}(kA1+1:end,:)]; 
            
            LC = [C_HSS{BOX.LAMBDA,ison}(:,1:kA1)*LA C_HSS{BOX.LAMBDA,ison}(:,kA1+1:end)*LB]; 
        end
        
        [Tout,Jout] = ID(RC,acc);
        [Tin,Jin] = ID(LC.',acc); 
        
        % Skeleton size k
        k = max(size(Tout,1),size(Tin,1));
        
        if ~(k == size(Tout,1))
            [Tout,Jout] = ID(RC, k);
        elseif ~(k == size(Tin,1))
            [Tin, Jin ] = ID(LC.',k);
        end
        
        % Record the outgoing skeletons:
        C_HSS{BOX.I_SKUP,ibox} = indskel_out(Jout(1:k));
        C_HSS{BOX.T_UP,ibox} = Tout;
        C_HSS{BOX.J_UP,ibox} = Jout;
        
        % Record the ingoing skeletons
        C_HSS{BOX.I_SKDN,ibox} = indskel_in(Jin(1:k));
        C_HSS{BOX.T_DN,ibox} = Tin;
        C_HSS{BOX.J_DN,ibox} = Jin;
        
        C_HSS{BOX.NSKEL,ibox} = length(indskel_out);
        C_HSS{BOX.KSKEL,ibox} = k;
        
        % Record pieces of R and L that need to be merged on the next level
        % rho = R^C(:,J(1:k))
        C_HSS{BOX.RHO,ibox} = RC(:,Jout(1:k));
        C_HSS{BOX.LAMBDA,ibox} = LC(Jin(1:k),:); 
        end
    end
end