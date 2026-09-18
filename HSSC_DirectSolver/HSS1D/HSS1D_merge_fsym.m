function [EHSS,X_sort,J] = HSS1D_merge_fsym(E1,E2,X,lev,acc,opts)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        [EHSS,X_sort,J] = HSS1D_merge_fsym(E1,E2,X,lev,acc,opts)
%
%    DESCRIPTION:
%        This function merges two fully-symmetric HSS structures encoded in E1 and E2 into a 2 x 2
%        block diagonal HSS matrix. In the HSS2D inversion algorithm, these blocks often come from
%        the restriction of Schur complement matrices.
%
%    INPUT:
%        E1,E2 are (50 x n1,n2) cell arrays containing HSS info for matrices E1,E2
%        X,lev are tree information relevant if we re-sort merged tree points in a curve. X is the
%            set of points where E1 and E2 are defined. The parity of lev determines how X is
%            oriented.
%        acc (double) is the desired accuracy.
%        opts is a string specifying whether or not to reduce each tree by one level.
%            opts == 'concatenate', we simply concatenate both trees.
%            opts == 'align', we align leaf boxes from both trees, which we assume have the
%                same structure up to sorting leaves appropriately. This is often needed if these
%                correspond to close-to-touching curves, for which the corresponding off-diagonal
%                interactions are not HSS unless we re-order the indices along the curve.
%
%    OUTPUT:
%        EHSS is a {50 x nboxes} cell array of the merged HSS matrix with nboxes = n1+n2+1.
%        X_sort is a (N x 2) array of the sorted set of points where EHSS is defined
%        J is a (N x 1) permutation vector such that X_sort = X(J,:).
%


global BOX 

n1 = size(E1,2);
n2 = size(E2,2); 

% STEP 1: Reduce one level of the tree, simplifying boxes with one child.
% This prevents the resulting HSS tree of becoming to deep / having small
% nodes. 

sm = 4;      

dif = abs(n1-n2);
while dif>0
  if n1<n2
      E2 = HSS1D_reducelev_fsym(E2,acc); 
      n2 = size(E2,2); 
  else  
      E1 = HSS1D_reducelev_fsym(E1,acc); 
      n1 = size(E1,2);
  end
  
  dif = abs(n1-n2);
end

if (n1>3 && n2>3) 
    cond1 = true; cond2 = true; 
    
    while (cond1 == 1 || cond2==1) && (n1>3 && n2>3)
        n1 = size(E1,2);
        n2 = size(E2,2);
        
        lf_small1 = 0; lf_small2 = 0; 
        lf_lson1 = 0; lf_lson2 = 0; 
    
        for ibox=max(n1,n2):-1:2
            if ibox<=n1
                
            if E1{BOX.C1,ibox}<=0 && E1{BOX.C2,ibox}<=0 
                if E1{BOX.END2,ibox} <= sm
                    lf_small1 = lf_small1+1; 
                end
                P = E1{BOX.PARENT,ibox}; 
                if E1{BOX.C1,P}<=0 || E1{BOX.C2,P}<=0
                    lf_lson1 = lf_lson1+1; 
                end
            else
                if E1{BOX.NSKEL,ibox} <= sm
                    lf_small1 = lf_small1+1; 
                end
            end
            end
       
            if ibox<=n2
                if E2{BOX.C1,ibox}<=0 && E2{BOX.C2,ibox}<=0 
                    if E2{BOX.END2,ibox} <= sm
                        lf_small2 = lf_small2+1; 
                    end
                    P = E2{BOX.PARENT,ibox}; 
                    if E2{BOX.C1,P}<=0 || E2{BOX.C2,P}<=0
                        lf_lson2 = lf_lson2+1; 
                    end
                else
                    if E2{BOX.NSKEL,ibox} <= sm      
                        lf_small2 = lf_small2+1;    
                    end
                end
            end
        end
       
        cond1 = (lf_lson1>0 || lf_small1>0);
        cond2 = (lf_lson2>0 || lf_small2>0);
    
        if cond1 == 1 || cond2==1
            E1 = HSS1D_reducelev_fsym(E1,acc);  
            E2 = HSS1D_reducelev_fsym(E2,acc); 
        end 
    end
end

% Depending on opts, HSS1D merge produces an HSS form for: 
if strcmp(opts,'concatenate') 
    % Concatenates the two trees, producing the matrix [E1 0 ; 0 E2] in
    % compressed form. 
    [EHSS,X_sort,J] = LOCAL_HSS1D_concatenate(E1,E2,X);
elseif strcmp(opts,'align')
    % Aligns the two trees, sorting X such that for leaf nodes, D = [D1 0 ; 0 D2]
    % L = [L1 0 ; 0 L2] and R = [R1 0 ; 0 R2]. 
    [EHSS,X_sort,J] = LOCAL_HSS1D_align(E1,E2,X,lev);
end

return

function [E,X_sort,J] = LOCAL_HSS1D_concatenate(E1,E2,X)
global BOX 

n1 = size(E1,2);
n2 = size(E2,2); 
nboxes = n1+n2+1; 
E = cell(size(E1,1),nboxes); 
ntot = E1{BOX.END2,1} + E2{BOX.END2,1}; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create root node: 
% Construct the top node.
E{BOX.LEVEL,1} = 0;
E{BOX.PARENT,1} = NaN;
E{BOX.C1,1} = 2;
E{BOX.C2,1} = 3;
E{BOX.END1,1} = 1;
E{BOX.END2,1} = ntot;

% roots of E1 and E2 are now children of new root node
E1{BOX.PARENT,1} = 0; 
E2{BOX.PARENT,1} = 0; 

depth = max(E1{BOX.LEVEL,n1},E2{BOX.LEVEL,n2}); 
numlev1 = zeros(1,depth+1); 
numlev2 = zeros(1,depth+1); 

for i=1:n1
   numlev1(E1{BOX.LEVEL,i}+1) = numlev1(E1{BOX.LEVEL,i}+1) + 1;  
   E1{BOX.LEVEL,i} = E1{BOX.LEVEL,i} + 1; 
end

for i=1:n2
   numlev2(E2{BOX.LEVEL,i}+1) = numlev2(E2{BOX.LEVEL,i}+1 ) + 1;
   E2{BOX.LEVEL,i} = E2{BOX.LEVEL,i} + 1; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Concatenate E = [E1 0 ; 0 E2] 

    n_j1 = 0; n_j2 = 0; n_jtot = 1; 
    for li = 1:depth+1
        if numlev1(li)>0
       % Modify info for E1 boxes
       i1_loc = (n_j1+1:n_j1+numlev1(li));
       
       for j = 1:numlev1(li)
          % Parent
          E1{BOX.PARENT,i1_loc(j)} = E1{BOX.PARENT,i1_loc(j)} + n_jtot; 
          
          if li>1
              E1{BOX.PARENT,i1_loc(j)} = E1{BOX.PARENT,i1_loc(j)} - n_j2 - numlev1(li-1); 
          end
          
          % Children (if the box is not a leaf)
          if li<E1{BOX.LEVEL,n1}
              if (E1{BOX.C1,i1_loc(j)}>0)
                 E1{BOX.C1,i1_loc(j)} = E1{BOX.C1,i1_loc(j)} + n_j2 + 1 + numlev2(li);
              end
              if (E1{BOX.C2,i1_loc(j)}>0)
                 E1{BOX.C2,i1_loc(j)} = E1{BOX.C2,i1_loc(j)} + n_j2 + 1 + numlev2(li); 
              end
          end
       end
       
       % Add E1 boxes to E on level li
       E(:,n_jtot+1:n_jtot+numlev1(li)) = E1(:,i1_loc); 
       n_j1   = n_j1  + numlev1(li); 
       n_jtot = n_jtot+ numlev1(li);
        end
       
       if numlev2(li)>0
       % Modify info for E2 boxes
       i2_loc = (n_j2+1:n_j2+numlev2(li)); 
       for j = 1:numlev2(li)
           % Parent
           E2{BOX.PARENT,i2_loc(j)} = E2{BOX.PARENT,i2_loc(j)} + n_jtot; 
           
           E2{BOX.PARENT,i2_loc(j)} = E2{BOX.PARENT,i2_loc(j)} - n_j1; 
          
           % Children (if the box is not a leaf)
           if li<E2{BOX.LEVEL,n2}
               if (E2{BOX.C1,i2_loc(j)}>0)
                  E2{BOX.C1,i2_loc(j)} = E2{BOX.C1,i2_loc(j)} + n_j1 + 1 + numlev1(li+1); 
               end
               if (E2{BOX.C2,i2_loc(j)}>0)
                  E2{BOX.C2,i2_loc(j)} = E2{BOX.C2,i2_loc(j)} + n_j1 + 1 + numlev1(li+1); 
               end
           end
           E2{BOX.END1,i2_loc(j)}  = E2{BOX.END1,i2_loc(j)}  + E1{BOX.END2,1}; % Update starting index
           E2{BOX.I_SKUP,i2_loc(j)} = E2{BOX.I_SKUP,i2_loc(j)} + E1{BOX.END2,1}; % indskel_out
           E2{BOX.I_SKDN,i2_loc(j)} = E2{BOX.I_SKDN,i2_loc(j)} + E1{BOX.END2,1}; % indskel_in
       end
       
       % Add E2 boxes to E on level li
       E(:,n_jtot+1:n_jtot+numlev2(li)) = E2(:,i2_loc); 
       n_j2   = n_j2  + numlev2(li); 
       n_jtot = n_jtot+ numlev2(li);
       end
       
    end
    
    % zero out level 1 interactions    
    for ibox = 3:(-1):2
        if size(E,2)>3
            ison1 = E{BOX.C1,ibox}; ison2 = E{BOX.C2,ibox};
            if ison1>0 && ison2>0
                E{BOX.I_SKUP,ibox} = [E{BOX.I_SKUP,ison1} E{BOX.I_SKUP,ison2}]; 
            else
                E{BOX.I_SKUP,ibox} = E{BOX.I_SKUP,max(ison1,ison2)};
            end 
        else
            E{BOX.I_SKUP,ibox} = E{BOX.END1,ibox} - 1 + (1:E{BOX.END2,ibox});
        end
        
        E{BOX.NSKEL,ibox} = length(E{BOX.I_SKUP,ibox});
        E{BOX.KSKEL,ibox} = 0;
        E{BOX.J_UP,ibox} = 1:E{BOX.NSKEL,ibox}; 
        E{BOX.T_UP,ibox} = zeros(0,E{BOX.NSKEL,ibox});
    end
    
    E{BOX.M_SIB,2} = zeros(0,0);
    E{BOX.M_SIB,3} = zeros(0,0);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X_sort = X; J = 1:size(X,1);  

return

function [E,X_sort,J] = LOCAL_HSS1D_align(E1,E2,X,lev)
global BOX 

n1 = size(E1,2);
nboxes = n1; 
E(1:5,:) = E1(1:5,:); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We assume the trees for E1 and E2 have the same structure, (they
% correspond to two sides of the same interface, or two curves close to 
% touching). Hence, we only need to match them and concatenate the 
% corresponding matrices. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sort and Match leaves (if necessary). This section can be commented out if
% both curves are already matched. 

Endopts1 = zeros(1,n1); Endopts2 = Endopts1;   
leaf = false(1,n1); 
% interface coordinate that is relevant depending on level
if mod(lev,2) == 1
    crd = 2; % vertical
else
    crd = 1; % horizontal
end

% Leaf box endpoints

for ibox = 1:n1
   if E{BOX.C1,ibox}<=0 & E{BOX.C2,ibox}<=0
       Endopts1(ibox) = X(E1{BOX.END1,ibox},crd); 
       Endopts2(ibox) = X(E2{BOX.END1,ibox}+E1{BOX.END2,1},crd); 
       leaf(ibox) = true; 
   end
end 

% Sort endpoints so that they match. 
Endopts1 = Endopts1(leaf); 
Endopts2 = Endopts2(leaf);    

J1_leaf = zeros(1,n1); J2_leaf = J1_leaf; 
[~,J1_leaf(leaf)] = sort(-Endopts1);
[~,J2_leaf(leaf)] = sort(-Endopts2);

J1_leaf(leaf) = 1:sum(leaf);        
%J2_leaf(leaf) = sum(leaf):-1:1;         
J2_leaf(leaf) = 1:sum(leaf);                                        
leafidx = 1:n1; leafidx = leafidx(leaf); 
ntot = E1{BOX.END2,1}+E2{BOX.END2,1}; 
J = zeros(1,ntot); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sorting of leaf nodes and recompression for non-leaf nodes

for ibox = nboxes:-1:1
    if (E{BOX.C1,ibox}<=0 & E{BOX.C2,ibox}<=0)
        ibox1 = leafidx(J1_leaf(ibox));
        ibox2 = leafidx(J2_leaf(ibox)); 
        
        E{BOX.END2,ibox} = E1{BOX.END2,ibox1} + E2{BOX.END2,ibox2}; 
        E{BOX.END1,ibox} = ntot - E{BOX.END2,ibox} + 1; 
        Ind = [(E1{BOX.END1,ibox1}-1+(1:E1{BOX.END2,ibox1})) ((E2{BOX.END1,ibox2}-1+(1:E2{BOX.END2,ibox2})) + E1{BOX.END2,1})];
                
        % Box Indices before subsampling
        indskel = E{BOX.END1,ibox}-1+(1:E{BOX.END2,ibox});
        J(indskel) = Ind;
        ntot = ntot - E{BOX.END2,ibox}; 
        
        E{BOX.NSKEL,ibox} = E{BOX.END2,ibox}; 
        
        % Interpolation Matrices
        m1 = E1{BOX.NSKEL,ibox1}; m2 = E2{BOX.NSKEL,ibox2}; 
        k1 = E1{BOX.KSKEL,ibox1}; k2 = E2{BOX.KSKEL,ibox2}; k = k1+k2; 
        
        R1 = [eye(k1)  E1{BOX.T_UP,ibox1}]; 
        R1(:,E1{BOX.J_UP,ibox1}) = R1;   
        R2 = [eye(k2)  E2{BOX.T_UP,ibox2}]; R2(:,E2{BOX.J_UP,ibox2}) = R2;
        R = [R1 zeros(k1,m2) ; zeros(k2,m1) R2]; 
        
        % Self Interactions                    
        E{BOX.M_SELF,ibox} = [E1{BOX.M_SELF,ibox1} zeros(m1,m2) ; zeros(m2,m1) E2{BOX.M_SELF,ibox2}]; 
        
        % Auxiliary info: 
        E{10,ibox} = ibox1; 
        E{11,ibox} = ibox2; 
    else
        % Box Indices before subsampling
        ison1 = E{BOX.C1,ibox};
        ison2 = E{BOX.C2,ibox};
        if (ison1>0 & ison2>0)
            indskel = [E{BOX.I_SKUP,ison1},E{BOX.I_SKUP,ison2}];
            E{BOX.END1,ibox} = min(E{BOX.END1,ison1},E{BOX.END1,ison2}); 
            E{BOX.END2,ibox} = E{BOX.END2,ison1} + E{BOX.END2,ison2};
            ibox1 = E1{BOX.PARENT,E{10,ison1}}; 
            ibox2 = E2{BOX.PARENT,E{11,ison1}}; 
        else
            ison = max([ison1 ison2]); 
            indskel = E{BOX.I_SKUP,ison};
            E{BOX.END1,ibox} = E{BOX.END1,ison}; 
            E{BOX.END2,ibox} = E{BOX.END2,ison}; 
            ibox1 = E1{BOX.PARENT,E{10,ison}}; 
            ibox2 = E2{BOX.PARENT,E{11,ison}}; 
        end
        
        %Factors from children recompression
        if (ison1>0 & ison2>0)
            R_c1 = E{BOX.RHO,ison1}; R_c2 = E{BOX.RHO,ison2};
            
            k1_1 = E1{BOX.KSKEL,E{10,ison1}}; k1_2 = E2{BOX.KSKEL,E{11,ison1}}; 
            k2_1 = E1{BOX.KSKEL,E{10,ison2}}; k2_2 = E2{BOX.KSKEL,E{11,ison2}};
        
            if ibox>1
            E{BOX.NSKEL,ibox} = length(indskel); 
            
            % Interpolation Matrices
            k1 = E1{BOX.KSKEL,ibox1}; k2 = E2{BOX.KSKEL,ibox2}; k = k1+k2; 
            R1 = [eye(k1)  E1{BOX.T_UP,ibox1}]; R1(:,E1{BOX.J_UP,ibox1}) = R1; 
            R2 = [eye(k2)  E2{BOX.T_UP,ibox2}]; R2(:,E2{BOX.J_UP,ibox2}) = R2;
            
            if E{10,ison1}<E{10,ison2}
                R11 = R1(:,1:k1_1); 
                R12 = R1(:,k1_1+1:end); 
            else
                R11 = R1(:,k2_1+1:end); 
                R12 = R1(:,1:k2_1); 
            end
            
            if E{11,ison1}<E{11,ison2}
                R21 = R2(:,1:k1_2); 
                R22 = R2(:,k1_2+1:end); 
            else
                R21 = R2(:,k2_2+1:end); 
                R22 = R2(:,1:k2_2); 
            end
            
            R = [R11*R_c1(1:k1_1,:) R12*R_c2(1:k2_1,:) ; R21*R_c1(k1_1+1:end,:) R22*R_c2(k2_1+1:end,:)]; 

            end
            
        else
            ison = max([ison1,ison2]); 
            E{BOX.NSKEL,ibox} = length(indskel); 
            
            % Interpolation Matrices
            k1 = E1{BOX.KSKEL,ibox1}; k2 = E2{BOX.KSKEL,ibox2}; k = k1+k2; 
            R1 = [eye(k1)  E1{BOX.T_UP,ibox1}]; R1(:,E1{BOX.J_UP,ibox1}) = R1; 
            R2 = [eye(k2)  E2{BOX.T_UP,ibox2}]; R2(:,E2{BOX.J_UP,ibox2}) = R2;
            
            k1_1 = E1{BOX.KSKEL,max(E1{4,ibox1},E1{BOX.C2,ibox1})};
            R_c = E{BOX.RHO,ison}; 
            R = [R1*R_c(1:k1_1,:) ; R2*R_c(k1_1+1:end,:)];   
         end
        
        % Sibling Interactions
        if (ison1>0 & ison2>0)
            B1 = E1{BOX.M_SIB,E{10,ison1}}; B2 = E2{BOX.M_SIB,E{11,ison1}};  
            
            %E{BOX.M_SIB,ison1} = [R_c1' R_c2']*[B1 zeros(size(B1,1),size(B2,2)) ; ...
            %              zeros(size(B2,1),size(B1,2)) B2]*[R_c1 ; R_c2];
            E{BOX.M_SIB,ison1} = R_c1(1:k1_1,:).'*B1*R_c2(1:k2_1,:) + R_c1(k1_1+1:end,:).'*B2*R_c2(k2_1+1:end,:);
            E{BOX.M_SIB,ison2} = E{BOX.M_SIB,ison1}.'; 
            
            %Clean Auxiliary info for sons
            E{BOX.RHO,ison1} = []; E{BOX.RHO,ison2} = []; 
            E{10,ison1} = []; E{11,ison1} = []; 
            E{10,ison2} = []; E{11,ison2} = [];
            clear R_c1 R_c2 R11 R12 R21 R22;
        else
            %Clean Auxiliary info for only son
            E{BOX.RHO,ison} = []; 
            E{10,ison} = []; E{11,ison} = [];
            clear R_c;
        end

        % Auxiliary info: 
        if ibox>1
            E{10,ibox} = ibox1; 
            E{11,ibox} = ibox2;
        end
    end
    
    if ibox>1
        % Recompression step (taking k = k1+k2, since R is block diagonal)
        [T,Jsk] = ID(R,k);
        k = size(T,1); 
        E{BOX.KSKEL,ibox} = k; 
        E{BOX.J_UP,ibox} = Jsk; %J
        E{BOX.I_SKUP,ibox} = indskel(Jsk(1:k)); %I_sk
        E{BOX.T_UP,ibox} = T; % T
    
        E{BOX.RHO,ibox} = R(:,Jsk(1:k)); 
    end
end
 
X_sort = X(J,:); 

return