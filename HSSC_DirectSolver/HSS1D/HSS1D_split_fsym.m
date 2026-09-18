function [K1,K2] = HSS1D_split_fsym(KHSS,I1,acc)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        [K1,K2] = HSS1D_split_fsym(KHSS,I1,acc)
%
%    DESCRIPTION:
%        This function splits a full-symmetric HSS matrix K into its diagonal blocks K1 = K(I1,I1)
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
%        K1,K2 are the HSS structures of matrices K1 = K(I1,I1), K2 = K(~I1,~I1)
%


global BOX 

nboxes = size(KHSS,2);
K1 = KHSS; 
K2 = K1; 
I1_src = cell(nboxes,1);
I1_sk = I1_src;
I2_src = cell(nboxes,1);
I2_sk = I2_src;
N1 = 0; N2 = 0;

if size(I1,1)>1
    I1 = I1'; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Count number of points on each leaf box and compute new endpoints
postorder = HSS1D_postorder(KHSS); 
for i = 1:nboxes
    ibox = postorder(i); 
    ison1 = KHSS{BOX.C1,ibox};
    ison2 = KHSS{BOX.C2,ibox};
    if ( ison1<=0 & ison2<=0 )
       I1_src{ibox} = I1((KHSS{BOX.END1,ibox} - 1 + (1:KHSS{BOX.END2,ibox})));
       I2_src{ibox} = ~I1_src{ibox}; 
       
       N1box = sum(I1_src{ibox}); N2box = sum(I2_src{ibox}); 
       K1{BOX.END1,ibox} = N1+1;       K2{BOX.END1,ibox} = N2+1;      
       K1{BOX.END2,ibox} = N1box;      K2{BOX.END2,ibox} = N2box;  
       N1 = N1+N1box;           N2 = N2+N2box; 
    elseif (ison1>0 & ison2>0)
       K1{BOX.END1,ibox} = min([K1{BOX.END1,ison1} K1{BOX.END1,ison2}]);
       K1{BOX.END2,ibox} = K1{BOX.END2,ison1}+K1{BOX.END2,ison2};
       K2{BOX.END1,ibox} = min([K2{BOX.END1,ison1} K2{BOX.END1,ison2}]);
       K2{BOX.END2,ibox} = K2{BOX.END2,ison1}+K2{BOX.END2,ison2};
    else
        ison = max([ison1 ison2]); 
        K1{BOX.END1,ibox} = K1{BOX.END1,ison};
        K1{BOX.END2,ibox} = K1{BOX.END2,ison};
        K2{BOX.END1,ibox} = K2{BOX.END1,ison};
        K2{BOX.END2,ibox} = K2{BOX.END2,ison};
    end
end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Step I: SPLIT and RECOMPRESS the two subtrees
for ibox = nboxes:(-1):1
    if ( (KHSS{BOX.C1,ibox}<=0) & (KHSS{BOX.C2,ibox}<=0) ) % ibox has no sons.
       %ind            = KHSS{BOX.END1,ibox} - 1 + (1:KHSS{BOX.END2,ibox});       
       ind1skel = K1{BOX.END1,ibox} - 1 + (1:K1{BOX.END2,ibox});
       ind2skel = K2{BOX.END1,ibox} - 1 + (1:K2{BOX.END2,ibox});
    else
        ison1       = KHSS{BOX.C1,ibox};
        ison2       = KHSS{BOX.C2,ibox};
        if (ison1>0 & ison2>0)  
            %I1_src{ibox} = [I1_sk{ison1} I1_sk{ison2}];
            %I2_src{ibox} = [I2_sk{ison1} I2_sk{ison2}];
           
            ind1skel = [K1{BOX.I_SKUP,ison1} K1{BOX.I_SKUP,ison2}];
            ind2skel = [K2{BOX.I_SKUP,ison1} K2{BOX.I_SKUP,ison2}];
        else
            ison = max([ison1 ison2]); 
            %I1_src{ibox} = I1_sk{ison};
            %I2_src{ibox} = I2_sk{ison};   
        
            ind1skel = K1{BOX.I_SKUP,ison};
            ind2skel = K2{BOX.I_SKUP,ison};
        end
    end
 
    if ibox>1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ison1 = KHSS{BOX.C1,ibox};
        ison2 = KHSS{BOX.C2,ibox};
        
        %Restricted Info for K1
        J = KHSS{BOX.J_UP,ibox}; k = KHSS{BOX.KSKEL,ibox}; 
        J_sk = J(1:k); 
        
        
        %if ison1>0 & ison2>0 
        %    K1{BOX.NSKEL,ibox} = K1{BOX.KSKEL,ison1} + K1{BOX.KSKEL,ison2}; %sum(I1_src{ibox}); % n_skel
        %elseif ison1>0 | ison2>0
        %    K1{BOX.NSKEL,ibox} = K1{BOX.KSKEL,max(ison1,ison2)};
        %else
        %    K1{BOX.NSKEL,ibox} = sum(I1_src{ibox}); 
        %end
        
        %K1{BOX.NSKEL,ibox} = sum(I1_src{ibox});          
        
        if (K1{BOX.END2,ibox}>0 & K1{BOX.END2,ibox}<KHSS{BOX.END2,ibox})  
            
            R = [eye(k)  KHSS{BOX.T_UP,ibox}];                   
                        
            if ( ison1<=0 & ison2<=0 )
                R(:,J) = R;
                R = R(:,I1_src{ibox}); 
            else
                R(:,J) = R;
                
                if (ison1>0 & ison2>0)
                    R = R*[K1{BOX.RHO,ison1} zeros(KHSS{BOX.KSKEL,ison1},K1{BOX.KSKEL,ison2}) ; zeros(KHSS{BOX.KSKEL,ison2},K1{BOX.KSKEL,ison1}) K1{BOX.RHO,ison2}];
                else
                    R = R*K1{BOX.RHO,ison}; 
                end
            end
            
            [T1,J1] = ID(R,acc);
            
            % Skeleton size k
            k1 = size(T1,1);
            
            K1{BOX.NSKEL,ibox} = length(J1); % n_skel 
            K1{BOX.KSKEL,ibox} = k1; %k_skel
            K1{BOX.J_UP,ibox} = J1; %J_up
            K1{BOX.I_SKUP,ibox} = ind1skel(J1(1:k1)); %I_skup
            K1{BOX.T_UP,ibox} = T1; % T_up(I_skup,I_rsup)          
            K1{BOX.RHO,ibox} = R(:,J1(1:k1)); 
            
        else
            if K1{BOX.END2,ibox} > 0
                I1_sk{ibox} = true(size(J_sk)); 
                K1{BOX.I_SKUP,ibox} = ind1skel(KHSS{BOX.J_UP,ibox}(1:KHSS{BOX.KSKEL,ibox}));
                K1{BOX.RHO,ibox} = eye(K1{BOX.KSKEL,ibox});
            else
                K1{BOX.NSKEL,ibox} = 0; 
                K1{BOX.KSKEL,ibox} = 0; 
                K1{BOX.RHO,ibox} = zeros(0,0); 
                %K1{BOX.J_UP,ibox} = zeros(1,0);   
                %K1{BOX.I_SKUP,ibox} = zeros(1,0); 
                %K1{BOX.T_UP,ibox} = zeros(0,0); 
                %}
            end
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Restricted Info for K2
        
        if (K2{BOX.END2,ibox}>0 & K2{BOX.END2,ibox}<KHSS{BOX.END2,ibox})   
            
            R = [eye(k)  KHSS{BOX.T_UP,ibox}];       
               
            if ( ison1<=0 & ison2<=0 )      
                R(:,J) = R;
                R = R(:,I2_src{ibox}); 
            else
                R(:,J) = R;
                
                if (ison1>0 & ison2>0)
                    R = R*[K2{BOX.RHO,ison1} zeros(KHSS{BOX.KSKEL,ison1},K2{BOX.KSKEL,ison2}) ; zeros(KHSS{BOX.KSKEL,ison2},K2{BOX.KSKEL,ison1}) K2{BOX.RHO,ison2}];
                else
                    R = R*K2{BOX.RHO,ison};
                end
            end
            
            [T2,J2] = ID(R,acc);
            
            % Skeleton size k
            k2 = size(T2,1);
        
            K2{BOX.NSKEL,ibox} = length(J2); % n_skel  
            K2{BOX.KSKEL,ibox} = k2; %k_skel
            K2{BOX.J_UP,ibox} = J2; %J_up
            K2{BOX.I_SKUP,ibox} = ind2skel(J2(1:k2)); %I_skup
            K2{BOX.T_UP,ibox} = T2; % T_up(I_skup,I_rsup)
            K2{BOX.RHO,ibox} = R(:,J2(1:k2)); 
            
        else
            if K2{BOX.END2,ibox} > 0
                I2_sk{ibox} = true(size(J_sk)); 
                K2{BOX.I_SKUP,ibox} = ind2skel(KHSS{BOX.J_UP,ibox}(1:KHSS{BOX.KSKEL,ibox}));
                K2{BOX.RHO,ibox} = eye(K2{BOX.KSKEL,ibox});
            else
                K2{BOX.NSKEL,ibox} = 0; 
                K2{BOX.KSKEL,ibox} = 0;
                K2{BOX.RHO,ibox} = zeros(0,0); 
            end
        end
    end
    
    % SELF and SIBLING INFO
    if ( (KHSS{BOX.C1,ibox}<=0) & (KHSS{BOX.C2,ibox}<=0) )
        %For leaf boxes, we restrict self interaction matrices
        K1{BOX.M_SELF,ibox} = KHSS{BOX.M_SELF,ibox}(I1_src{ibox},I1_src{ibox}); 
        K2{BOX.M_SELF,ibox} = KHSS{BOX.M_SELF,ibox}(I2_src{ibox},I2_src{ibox});
    
    else
        %For non-leaf, we restrict sibling interactions
        ison1 = KHSS{BOX.C1,ibox};
        ison2 = KHSS{BOX.C2,ibox};     
        
        if (ison1>0 & ison2>0)
        %if (K1{BOX.NSKEL,ison1}>0 & K1{BOX.NSKEL,ison1}<KHSS{BOX.NSKEL,ison1}) & (K1{BOX.NSKEL,ison2}>0 & K1{BOX.NSKEL,ison2}<KHSS{BOX.NSKEL,ison2})
            if (K1{BOX.END2,ison1}>0 ) & (K1{BOX.END2,ison2}>0)
                K1{BOX.M_SIB,ison1} = K1{BOX.RHO,ison1}.'*KHSS{BOX.M_SIB,ison1}*K1{BOX.RHO,ison2}; 
                K1{BOX.M_SIB,ison2} = K1{BOX.M_SIB,ison1}.'; 
            end
       
        %if (K2{BOX.NSKEL,ison1}>0 & K2{BOX.NSKEL,ison1}<KHSS{BOX.NSKEL,ison1}) & (K2{BOX.NSKEL,ison2}>0 & K2{BOX.NSKEL,ison2}<KHSS{BOX.NSKEL,ison2})
            if (K2{BOX.END2,ison1}>0) & (K2{BOX.END2,ison2}>0)
                K2{BOX.M_SIB,ison1} = K2{BOX.RHO,ison1}.'*KHSS{BOX.M_SIB,ison1}*K2{BOX.RHO,ison2}; 
                K2{BOX.M_SIB,ison2} = K2{BOX.M_SIB,ison1}.';  
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
    if K1{BOX.END2,j} > 0
        full1(j) = 1;
        sum_full1(j) = sum(full1(1:j)); 
    end
    if K2{BOX.END2,j} > 0   
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
while ((nc1==0 || nc2==0) & ncsq>0)
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
while ((nc1==0 || nc2==0) & ncsq>0)
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
            K1{BOX.PARENT,K1{BOX.C1,j}} = j; 
        else
            K1{BOX.C1,j} = -1; 
        end
    end
    if K1{BOX.C2,j}>0
        K1{BOX.C2,j} = sum_full1(K1{BOX.C2,j});
        if K1{BOX.C2,j}>0
            K1{BOX.PARENT,K1{BOX.C2,j}} = j; 
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
            K2{BOX.PARENT,K2{BOX.C1,j}} = j; 
        else
            K2{BOX.C1,j} = -1; 
        end
    end
    if K2{BOX.C2,j}>0
        K2{BOX.C2,j} = sum_full2(K2{BOX.C2,j});
        if K2{BOX.C2,j}>0
            K2{BOX.PARENT,K2{BOX.C2,j}} = j; 
        else
            K2{BOX.C2,j} = -1;
        end
    end
end
       
return