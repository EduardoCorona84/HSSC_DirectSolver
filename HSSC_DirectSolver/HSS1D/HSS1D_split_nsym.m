function [K1,K2] = HSS1D_split_nsym(KHSS,I1_dn,I1_up,acc)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        [K1,K2] = HSS1D_split_nsym(KHSS,I1_dn,I1_up,acc)
%
%    DESCRIPTION:
%        This function splits a non-symmetric HSS matrix K into its diagonal blocks K1 = K(I1,I1)
%        and K2 = K(~I1,~I1). In order to do this, we note that:
%                K(I,I) = D_n(I,I) + R(:,I)'[K^{n-1}]R(:,I)
%        So, leaf nodes must be restricted and the rest re-compressed. If all points are removed
%        (zero) from a certain subtree, the whole subtree is removed. If the root node has only one
%        child it is removed. This process is repeated as many times as necessary.
%        This function is used for constructing auxiliary matrices E,F^-1 in HSS2D inversion.
%
%    INPUT:
%        KHSS is the HSS structure of matrix K.
%        I1 is a (1 x N) boolean vector specifying which points we pick.
%        acc (double) is the desired accuracy.
%
%    OUTPUT:
%        K1,K2 is the HSS structure of matrices K(I1,I1), K(~I1,~I1).
%


global BOX 

nboxes = size(KHSS,2);
K1 = KHSS; 
K2 = K1; 
I1_srcup = cell(nboxes,1);
I1_srcdn = I1_srcup; 
I1_skup = I1_srcup;
I1_skdn = I1_srcup; 
I2_srcup = cell(nboxes,1);
I2_srcdn = I2_srcup; 
I2_skup = I2_srcup;
I2_skdn = I2_srcup;
N1 = 0; N2 = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Count number of points on each leaf box and compute new endpoints
for ibox = 1:nboxes
    if ( (KHSS{BOX.C1,ibox}<=0) && (KHSS{BOX.C2,ibox}<=0) )
       I1_srcup{ibox} = I1_up((KHSS{BOX.END1,ibox} - 1 + (1:KHSS{BOX.END2,ibox})));
       I1_srcdn{ibox} = I1_dn((KHSS{BOX.END1,ibox} - 1 + (1:KHSS{BOX.END2,ibox}))); 
       I2_srcup{ibox} = ~I1_srcup{ibox}; 
       I2_srcdn{ibox} = ~I1_srcup{ibox}; 
       
       N1box = sum(I1_srcup{ibox}); N2box = sum(~I1_srcup{ibox}); 
       K1{BOX.END1,ibox} = N1+1;       K2{BOX.END1,ibox} = N2+1;      
       K1{BOX.END2,ibox} = N1box;      K2{BOX.END2,ibox} = N2box;  
       N1 = N1+N1box;           N2 = N2+N2box; 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Step I: SPLIT and RECOMPRESS the two subtrees
for ibox = nboxes:(-1):1
    if ( (KHSS{BOX.C1,ibox}<=0) && (KHSS{BOX.C2,ibox}<=0) ) % ibox has no sons.
       %ind            = KHSS{BOX.END1,ibox} - 1 + (1:KHSS{BOX.END2,ibox});       
       ind1 = K1{BOX.END1,ibox} - 1 + (1:K1{BOX.END2,ibox});
       ind1skel_out = ind1; ind1skel_in = ind1; 
       ind2 = K2{BOX.END1,ibox} - 1 + (1:K2{BOX.END2,ibox});
       ind2skel_out = ind2; ind2skel_in = ind2;
    else
        ison1       = KHSS{BOX.C1,ibox};
        ison2       = KHSS{BOX.C2,ibox};
        if (ison1>0 && ison2>0)
            I1_srcup{ibox} = [I1_skup{ison1},I1_skup{ison2}];
            I1_srcdn{ibox} = [I1_skdn{ison1},I1_skdn{ison2}];
            I2_srcup{ibox} = [I2_skup{ison1},I2_skup{ison2}];
            I2_srcdn{ibox} = [I2_skdn{ison1},I2_skdn{ison2}];
            
            K1{BOX.END1,ibox} = K1{BOX.END1,ison1}; K1{BOX.END2,ibox} = K1{BOX.END2,ison1}+K1{BOX.END2,ison2}; 
            K2{BOX.END1,ibox} = K2{BOX.END1,ison1}; K2{BOX.END2,ibox} = K2{BOX.END2,ison1}+K2{BOX.END2,ison2};
        
            ind1skel_out = [K1{BOX.I_SKUP,ison1},K1{BOX.I_SKUP,ison2}];
            ind1skel_in  = [K1{BOX.I_SKDN,ison1},K1{BOX.I_SKDN,ison2}];
            ind2skel_out = [K2{BOX.I_SKUP,ison1},K2{BOX.I_SKUP,ison2}];
            ind2skel_in  = [K2{BOX.I_SKDN,ison1},K2{BOX.I_SKDN,ison2}];
        else
            ison = max([ison1 ison2]); 
            I1_srcup{ibox} = I1_skup{ison};
            I1_srcdn{ibox} = I1_skdn{ison};
            I2_srcup{ibox} = I2_skup{ison};
            I2_srcdn{ibox} = I2_skdn{ison};
            
            K1{BOX.END1,ibox} = K1{BOX.END1,ison}; K1{BOX.END2,ibox} = K1{BOX.END2,ison}; 
            K2{BOX.END1,ibox} = K2{BOX.END1,ison}; K2{BOX.END2,ibox} = K2{BOX.END2,ison};
        
            ind1skel_out = K1{BOX.I_SKUP,ison};
            ind1skel_in  = K1{BOX.I_SKDN,ison};
            ind2skel_out = K2{BOX.I_SKUP,ison};
            ind2skel_in  = K2{BOX.I_SKDN,ison};
        end
        
        
      
    end
 
    if ibox>1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Restricted Info for K1
        J_up = KHSS{BOX.J_UP,ibox}; k = KHSS{BOX.KSKEL,ibox}; 
        J_skup = J_up(1:k); 
        J_dn = KHSS{BOX.J_DN,ibox}; J_skdn = J_dn(1:k);
        K1{BOX.NSKEL,ibox} = sum(I1_srcup{ibox}); % n_skel
        if (K1{BOX.NSKEL,ibox}>0 && K1{BOX.NSKEL,ibox}<KHSS{BOX.NSKEL,ibox})
            
            R = [eye(k)  KHSS{BOX.T_UP,ibox}];       
            L = [eye(k) ; KHSS{BOX.T_DN,ibox}.']; 
            
                        
            if ( (KHSS{BOX.C1,ibox}<=0) && (KHSS{BOX.C2,ibox}<=0) )
                R(:,J_up) = R;
                R = R(:,I1_srcup{ibox}); 
                
                L(J_dn,:) = L; 
                L = L(I1_srcdn{ibox},:); 
            else
                R(:,J_up) = R;
                L(J_dn,:) = L; 
                ison1 = KHSS{BOX.C1,ibox};
                ison2 = KHSS{BOX.C2,ibox};
                
                if (ison1>0 && ison2>0)
                    R = R*[K1{BOX.RHO,ison1} zeros(KHSS{BOX.KSKEL,ison1},K1{BOX.KSKEL,ison2}) ; zeros(KHSS{BOX.KSKEL,ison2},K1{BOX.KSKEL,ison1}) K1{BOX.RHO,ison2}];
                    L = [K1{BOX.LAMBDA,ison1} zeros(K1{BOX.KSKEL,ison1},KHSS{BOX.KSKEL,ison2}); zeros(K1{BOX.KSKEL,ison2},KHSS{BOX.KSKEL,ison1}) K1{BOX.LAMBDA,ison2}]*L; 
                else
                    R = R*K1{BOX.RHO,ison}; 
                    L = K1{BOX.LAMBDA,ison}*L; 
                end
            end
            
            [T1_up,J1_up] = ID(R,acc);
            [T1_dn,J1_dn] = ID(L.',acc);
            
            % Skeleton size k
            k1 = max(size(T1_up,1),size(T1_dn,1));
        
            % We sometimes need to enforce that the outgoing rank = incoming rank.
            % (Note that this is implemented rather clumsily.)
            if ~(k1 == size(T1_up,1))
                [T1_up,J1_up] = ID(R, k1);
            elseif ~(k1 == size(T1_dn,1))
                [T1_dn, J1_dn ] = ID(L.',k1);
            end
            
            K1{BOX.NSKEL,ibox} = length(J1_up); % n_skel
            K1{BOX.KSKEL,ibox} = k1; %k_skel
            K1{BOX.J_UP,ibox} = J1_up; %J_up
            K1{BOX.I_SKUP,ibox} = ind1skel_out(J1_up(1:k1)); %I_skup
            K1{BOX.T_UP,ibox} = T1_up; % T_up(I_skup,I_rsup)
            
            I1_skup{ibox} = I1_srcup{ibox}(J_skup); 
            
            K1{BOX.J_DN,ibox} = J1_dn; %J_dn
            K1{BOX.I_SKDN,ibox} = ind1skel_in(J1_dn(1:k1)); %I_skdn
            K1{BOX.T_DN,ibox} = T1_dn; % T_dn(I_skdn,I_rsdn)
            
            I1_skdn{ibox} = I1_srcdn{ibox}(J_skdn); 
            
            K1{BOX.RHO,ibox} = R(:,J1_up(1:k1)); 
            K1{BOX.LAMBDA,ibox} = L(J1_dn(1:k1),:); 
            
        else
            if K1{BOX.NSKEL,ibox} > 0
                I1_skup{ibox} = true(size(J_skup)); 
                I1_skdn{ibox} = true(size(J_skdn)); 
                K1{BOX.I_SKUP,ibox} = ind1skel_out(KHSS{BOX.J_UP,ibox}(1:KHSS{BOX.KSKEL,ibox}));
                K1{BOX.I_SKDN,ibox} = ind1skel_in(KHSS{BOX.J_DN,ibox}(1:KHSS{BOX.KSKEL,ibox})); 
                K1{BOX.RHO,ibox} = eye(K1{BOX.KSKEL,ibox}); 
                K1{BOX.LAMBDA,ibox} = K1{BOX.RHO,ibox}; 
            else
                K1{BOX.KSKEL,ibox} = 0; 
                I1_skup{ibox} = false(size(J_skup)); 
                I1_skdn{ibox} = false(size(J_skdn)); 
                K1{BOX.RHO,ibox} = sparse(K1{BOX.KSKEL,ibox},0); 
                K1{BOX.LAMBDA,ibox} = sparse(0,K1{BOX.KSKEL,ibox}); 
            end
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Restricted Info for K2
        
        K2{BOX.NSKEL,ibox} = sum(I2_srcup{ibox}); % n_skel
        if (K2{BOX.NSKEL,ibox}>0 && K2{BOX.NSKEL,ibox}<KHSS{BOX.NSKEL,ibox})
            
            R = [eye(k)  KHSS{BOX.T_UP,ibox}];       
            L = [eye(k) ; KHSS{BOX.T_DN,ibox}.']; 
            
                        
            if ( (KHSS{BOX.C1,ibox}<=0) && (KHSS{BOX.C2,ibox}<=0) )
                R(:,J_up) = R;
                R = R(:,I2_srcup{ibox}); 
                
                L(J_dn,:) = L; 
                L = L(I2_srcdn{ibox},:); 
            else
                R(:,J_up) = R;
                L(J_dn,:) = L; 
                ison1 = KHSS{BOX.C1,ibox};
                ison2 = KHSS{BOX.C2,ibox};
                
                if (ison1>0 && ison2>0)
                    R = R*[K2{BOX.RHO,ison1} zeros(KHSS{BOX.KSKEL,ison1},K2{BOX.KSKEL,ison2}) ; zeros(KHSS{BOX.KSKEL,ison2},K2{BOX.KSKEL,ison1}) K2{BOX.RHO,ison2}];
                    L = [K2{BOX.LAMBDA,ison1} zeros(K2{BOX.KSKEL,ison1},KHSS{BOX.KSKEL,ison2}); zeros(K2{BOX.KSKEL,ison2},KHSS{BOX.KSKEL,ison1}) K2{BOX.LAMBDA,ison2}]*L; 
                else
                    R = R*K2{BOX.RHO,ison}; 
                    L = K2{BOX.LAMBDA,ison}*L; 
                end
            end
            
            [T2_up,J2_up] = ID(R,acc);
            [T2_dn,J2_dn] = ID(L.',acc);
            
            % Skeleton size k
            k2 = max(size(T2_up,1),size(T2_dn,1));
        
            % We sometimes need to enforce that the outgoing rank = incoming rank.
            % (Note that this is implemented rather clumsily.)
            if ~(k2 == size(T2_up,1))
                [T2_up,J2_up] = ID(R, k2);
            elseif ~(k2 == size(T2_dn,1))
                [T2_dn, J2_dn ] = ID(L.',k2);
            end
            
            K2{BOX.NSKEL,ibox} = length(J2_up); % n_skel
            K2{BOX.KSKEL,ibox} = k2; %k_skel
            K2{BOX.J_UP,ibox} = J2_up; %J_up
            K2{BOX.I_SKUP,ibox} = ind2skel_out(J2_up(1:k2)); %I_skup
            K2{BOX.T_UP,ibox} = T2_up; % T_up(I_skup,I_rsup)
            
            I2_skup{ibox} = I2_srcup{ibox}(J_skup); 
            
            K2{BOX.J_DN,ibox} = J2_dn; %J_dn
            K2{BOX.I_SKDN,ibox} = ind2skel_in(J2_dn(1:k2)); %I_skdn
            K2{BOX.T_DN,ibox} = T2_dn; % T_dn(I_skdn,I_rsdn)
            
            I2_skdn{ibox} = I2_srcdn{ibox}(J_skdn); 
            
            K2{BOX.RHO,ibox} = R(:,J2_up(1:k2)); 
            K2{BOX.LAMBDA,ibox} = L(J2_dn(1:k2),:); 
            
        else
            if K2{BOX.NSKEL,ibox} > 0
                I2_skup{ibox} = true(size(J_skup)); 
                I2_skdn{ibox} = true(size(J_skdn));
                K2{BOX.I_SKUP,ibox} = ind2skel_out(KHSS{BOX.J_UP,ibox}(1:KHSS{BOX.KSKEL,ibox}));
                K2{BOX.I_SKDN,ibox} = ind2skel_in(KHSS{BOX.J_DN,ibox}(1:KHSS{BOX.KSKEL,ibox})); 
                K2{BOX.RHO,ibox} = eye(K2{BOX.KSKEL,ibox}); 
                K2{BOX.LAMBDA,ibox} = K2{BOX.RHO,ibox}; 
            else
                I2_skup{ibox} = false(size(J_skup)); 
                I2_skdn{ibox} = false(size(J_skdn)); 
                K2{BOX.KSKEL,ibox} = 0;
                K2{BOX.RHO,ibox} = sparse(K2{BOX.KSKEL,ibox},0); 
                K2{BOX.LAMBDA,ibox} = sparse(0,K2{BOX.KSKEL,ibox}); 
            end
        end
    end
    
    % SELF and SIBLING INFO
    if ( (KHSS{BOX.C1,ibox}<=0) && (KHSS{BOX.C2,ibox}<=0) )
        %For leaf boxes, we restrict self interaction matrices
        K1{BOX.M_SELF,ibox} = KHSS{BOX.M_SELF,ibox}(I1_srcdn{ibox},I1_srcup{ibox}); 
        K2{BOX.M_SELF,ibox} = KHSS{BOX.M_SELF,ibox}(I2_srcdn{ibox},I2_srcup{ibox});
    
    else
        %For non-leaf, we restrict sibling interactions
        ison1 = KHSS{BOX.C1,ibox};
        ison2 = KHSS{BOX.C2,ibox};     
        
        if (ison1>0 && ison2>0)
        %if (K1{BOX.NSKEL,ison1}>0 && K1{BOX.NSKEL,ison1}<KHSS{BOX.NSKEL,ison1}) && (K1{BOX.NSKEL,ison2}>0 && K1{BOX.NSKEL,ison2}<KHSS{BOX.NSKEL,ison2})
            if (K1{BOX.NSKEL,ison1}>0 ) && (K1{BOX.NSKEL,ison2}>0)
                K1{BOX.M_SIB,ison1} = K1{BOX.LAMBDA,ison1}*KHSS{BOX.M_SIB,ison1}*K1{BOX.RHO,ison2}; 
                K1{BOX.M_SIB,ison2} = K1{BOX.LAMBDA,ison2}*KHSS{BOX.M_SIB,ison2}*K1{BOX.RHO,ison1}; 
            end
       
        %if (K2{BOX.NSKEL,ison1}>0 && K2{BOX.NSKEL,ison1}<KHSS{BOX.NSKEL,ison1}) && (K2{BOX.NSKEL,ison2}>0 && K2{BOX.NSKEL,ison2}<KHSS{BOX.NSKEL,ison2})
            if (K2{BOX.NSKEL,ison1}>0) && (K2{BOX.NSKEL,ison2}>0)
                K2{BOX.M_SIB,ison1} = K2{BOX.LAMBDA,ison1}*KHSS{BOX.M_SIB,ison1}*K2{BOX.RHO,ison2}; 
                K2{BOX.M_SIB,ison2} = K2{BOX.LAMBDA,ison2}*KHSS{BOX.M_SIB,ison2}*K2{BOX.RHO,ison1}; 
            end
        end
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLEAN UP of K1 and K2 (remove empty nodes and spurious root nodes)
full1 = zeros(1,nboxes); full2 = full1; 
sum_full1 = full1; sum_full2 = full2; 

% Check which boxes have points in them
full1(1) = 1; full2(1) = 1; 
sum_full1(1) = 1; sum_full2(1) = 1; 
for j=2:nboxes
    if K1{BOX.NSKEL,j} > 0
        full1(j) = 1;
        sum_full1(j) = sum(full1(1:j)); 
    end
    if K2{BOX.NSKEL,j} > 0
        full2(j) = 1;
        sum_full2(j) = sum(full2(1:j)); 
    end
end

% Root node(s) that need to be deleted for K1
mlev1 = 0; 
j = 1; 
nc1 = K1{BOX.NSKEL,K1{4,j}}; 
nc2 = K1{BOX.NSKEL,K1{5,j}};
ncsq = 1; 
while ((nc1==0 || nc2==0) && ncsq>0)
   mlev1 = mlev1+1; 
   full1(j) = 0; 
   if K1{BOX.NSKEL,K1{4,j}} > K1{BOX.NSKEL,K1{5,j}}
       j = K1{BOX.C1,j}; 
   else
       j = K1{BOX.C2,j}; 
   end
   
   if K1{BOX.C1,j}>0 
       nc1 = K1{BOX.NSKEL,K1{4,j}};
   else
       nc1 = 0; 
   end
   
   if K1{BOX.C2,j}>0 
       nc2 = K1{BOX.NSKEL,K1{5,j}};
   else
       nc2 = 0; 
   end
   
   ncsq = nc1^2+nc2^2; 
end

% Root node(s) that need to be deleted for K2
mlev2 = 0; 
j = 1; 
nc1 = K2{BOX.NSKEL,K2{4,j}}; 
nc2 = K2{BOX.NSKEL,K2{5,j}};
ncsq = 1; 
while ((nc1==0 || nc2==0) && ncsq>0)
   mlev2 = mlev2+1; 
   full2(j) = 0; 
   if K2{BOX.NSKEL,K2{4,j}} > K2{BOX.NSKEL,K2{5,j}}
       j = K2{BOX.C1,j}; 
   else
       j = K2{BOX.C2,j}; 
   end
   
   if K2{BOX.C1,j}>0 
       nc1 = K2{BOX.NSKEL,K2{4,j}};
   else
       nc1 = 0; 
   end
   
   if K2{BOX.C2,j}>0 
       nc2 = K2{BOX.NSKEL,K2{5,j}};
   else
       nc2 = 0; 
   end
   
   ncsq = nc1^2+nc2^2; 
end

sum_full1 = sum_full1 - mlev1*(sum_full1>0);
sum_full2 = sum_full2 - mlev2*(sum_full2>0);

% DELETE 'empty' rows and reassign parent and children information
K1 = K1(:,full1==1); 
K2 = K2(:,full2==1); 

K1{BOX.PARENT,1} = NaN; K2{BOX.PARENT,1} = NaN; 

for j=1:size(K1,2)
    % Reduce levels 
    K1{BOX.LEVEL,j} = K1{BOX.LEVEL,j} - mlev1; 
    %New Parent/Children Info
    if K1{BOX.C1,j}>0
        K1{BOX.C1,j} = sum_full1(K1{BOX.C1,j});
        if K1{BOX.C1,j}>0
            K1{BOX.PARENT,K1{4,j}} = j; 
        else
            K1{BOX.C1,j} = -1; 
        end
    end
    if K1{BOX.C2,j}>0
        K1{BOX.C2,j} = sum_full1(K1{BOX.C2,j});
        if K1{BOX.C2,j}>0
            K1{BOX.PARENT,K1{5,j}} = j; 
        else
            K1{BOX.C2,j} = -1; 
        end
    end
end


for j=1:size(K2,2)
    % Reduce levels 
    K2{BOX.LEVEL,j} = K2{BOX.LEVEL,j} - mlev2; 
    %New Parent/Children Info
    if K2{BOX.C1,j}>0
        K2{BOX.C1,j} = sum_full2(K2{BOX.C1,j});
        if K2{BOX.C1,j}>0
            K2{BOX.PARENT,K2{4,j}} = j; 
        else
            K2{BOX.C1,j} = -1; 
        end
    end
    if K2{BOX.C2,j}>0
        K2{BOX.C2,j} = sum_full2(K2{BOX.C2,j});
        if K2{BOX.C2,j}>0
            K2{BOX.PARENT,K2{5,j}} = j; 
        else
            K2{BOX.C2,j} = -1;
        end
    end
end
    
       
return