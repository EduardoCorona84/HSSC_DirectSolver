function A_NEW = HSS1D_recompress_fsym(A,acc)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        A_NEW = HSS1D_recompress_fsym(A,acc)
%
%    DESCRIPTION:
%        This function recompresses a full-symmetric HSS matrix A, obtaining an optimal
%        representation (given a goal accuracy) by recursively re-compressing factors that form
%        neutered columns and rows. This function should be used when A is suboptimal due to HSS
%        arithmetic (e.g. sum, restriction, low rank update)
%
%        (1) The first part of this algorithm uses SVDs to compress A, a method a adapted from
%            Jianlin Xia's algorithm on "On the complexity of somehierarchical structured matrix
%            algorithms", p. 14-19.
%        (2) We then perform an upward pass where we re-format back into HSS1D form; we define new
%            skeleton information and interpolation operators and update sibling interactions
%            accordingly.
%
%    INPUT:
%        A is a {50 x nboxes} cell array containing HSS info for matrix A. A{36,ibox} and
%            A{37,ibox} are used to store auxiliary matrices \tilde{T} and R_hat / \tilde{R} as
%            defined in the paper.
%        acc (double) is the desired accuracy
%
%    OUTPUT:
%        A_NEW is the {50 x nbox_NEW} cell array containing HSS info for recompressed A.
%


global BOX 

nboxes = size(A,2); 
A_NEW(BOX.DATA,:) = A(BOX.DATA,:); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (I) Downward pass, Compression of neutered columns / rows using SVDs
for ibox = 1:nboxes
    if ( (A{BOX.C1,ibox}>0) || (A{BOX.C2,ibox}>0) ) % ibox has at least one son
        ison1       = A{BOX.C1,ibox};
        ison2       = A{BOX.C2,ibox};
        if (ison1>0 && ison2>0)

            k1 = A{BOX.KSKEL,ison1}; R1 = [eye(k1) A{BOX.T_UP,ison1}]; R1(:,A{BOX.J_UP,ison1}) = R1;
            k2 = A{BOX.KSKEL,ison2}; R2 = [eye(k2) A{BOX.T_UP,ison2}]; R2(:,A{BOX.J_UP,ison2}) = R2;
            
            Tson1 = A{BOX.M_SIB,ison1}; 
            Tson2 = Tson1.';
            if ibox==1 
                [P,D,Q] = svd(Tson2,0);
                k1_NEW = sum(diag(D)>acc); 
                P = P(:,1:k1_NEW); D = D(1:k1_NEW,1:k1_NEW); Q = Q(:,1:k1_NEW); 
                
                % Sibling interactions 
                A_NEW{BOX.M_SIB,ison2} = D; A_NEW{BOX.M_SIB,ison1} = D; 
                % \tilde{T}_i
                A_NEW{36,ison1} = D; A_NEW{36,ison2} = D; 
                % \hat{R}_i 
                A_NEW{37,ison1} = Q'*R1; 
                A_NEW{37,ison2} = P.'*R2; 
            
            else
               % \tilde{T} for parent and R_hat
               Tp = A_NEW{36,ibox};   
               Rp = A_NEW{37,ibox};
               
               % Define \tilde{T} for children
               Tson1 = [Tson1 ; Tp*Rp(:,k1+1:end)];
               Tson2 = [Tson2 ; Tp*Rp(:,1:k1)]; 
               
               % clear Tp
               A_NEW{36,ibox} = [];
               
               % SVD(T_i)
                [Y1,D1,Q1] = svd(Tson2,0); k1_NEW = sum(diag(D1)>acc); 
                [Y2,D2,Q2] = svd(Tson1,0); k2_NEW = sum(diag(D2)>acc); 
            
                % Re-size the factors of the svd
                Y1 = Y1(:,1:k1_NEW); D1 = D1(1:k1_NEW,1:k1_NEW); Q1 = Q1(:,1:k1_NEW);
                Y2 = Y2(:,1:k2_NEW); D2 = D2(1:k2_NEW,1:k2_NEW); Q2 = Q2(:,1:k2_NEW);
            
                % Sibling Interactions
                A_NEW{BOX.M_SIB,ison1} = Q1.'*Y2(1:k1,:)*D2; 
                A_NEW{BOX.M_SIB,ison2} = Q2.'*Y1(1:k2,:)*D1; 
                % \tilde{T} (diagonal matrix)
                A_NEW{36,ison1} = D1; A_NEW{36,ison2} = D2; 
                % \hat{R}_i 
                A_NEW{37,ison1} = Q1'*R1; A_NEW{37,ison2} = Q2'*R2;
                % \tilde{R}_p (new "interpolation" operator)
                A_NEW{37,ibox} = [Rp(:,1:k1)*Q1 Rp(:,k1+1:end)*Q2]; 
            end
            
        else
            ison = max(ison1,ison2);
            k = A{BOX.KSKEL,ison}; R = [eye(k) A{BOX.T_UP,ison}]; R(:,A{BOX.J_UP,ison}) = R;
            Tp = A_NEW{36,ibox}; Rp = A_NEW{37,ibox};
            
            % Since there are no sibling interactions, this box interacts
            % with the neutered column/row only through its parent. 
            Tson = Tp*Rp; 
            
            % SVD(T_i)
            [~,D,Q] = svd(Tson,0); k_NEW = sum(diag(D)>acc); 
            % Re-size the factors of the svd
            D = D(1:k_NEW,1:k_NEW); Q = Q(:,1:k_NEW);
            % \tilde{T} 
            A_NEW{36,ison} = D;
            % \hat{R}_i 
            A_NEW{37,ison} = Q'*R;
            % \tilde{R}_p 
            A_NEW{37,ibox} = Rp*Q; 
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (II) Re-formatting to ID HSS: 

for ibox = nboxes:(-1):1
    if ( (A{BOX.C1,ibox}<=0) && (A{BOX.C2,ibox}<=0) ) % ibox has no sons.
        %indices
        indskel = A{BOX.END1,ibox} - 1 + (1:A{BOX.END2,ibox});
        % Self-interactions are the same
        A_NEW{BOX.M_SELF,ibox} = A{BOX.M_SELF,ibox};
        % We need to re-format the new interpolation operator
        RNEW = A_NEW{37,ibox}; 
    else
        ison1       = A_NEW{BOX.C1,ibox};
        ison2       = A_NEW{BOX.C2,ibox};
        if (ison1>0 && ison2>0)
            indskel = [A_NEW{BOX.I_SKUP,ison1},A_NEW{BOX.I_SKUP,ison2}];
            
             % Sibling Interaction Matrices
            A_NEW{BOX.M_SIB,ison1} = A_NEW{BOX.RHO,ison1}.'*A_NEW{BOX.M_SIB,ison1}*A_NEW{BOX.RHO,ison2};
            A_NEW{BOX.M_SIB,ison2} = A_NEW{BOX.M_SIB,ison1}.'; 
        else
            ison = max([ison1,ison2]); 
            indskel = A_NEW{BOX.I_SKUP,ison};
        end
        
        if ibox>1
        
            % Right / Outgoing Interpolation Matrices 
            RA = A_NEW{37,ibox}; 
            
            % We include factors coming from recompression of leaves. 
            if (ison1>0 && ison2>0)
                kA1 = size(A_NEW{BOX.RHO,ison1},1); 
                RNEW = [RA(:,1:kA1)*A_NEW{BOX.RHO,ison1} RA(:,kA1+1:end)*A_NEW{BOX.RHO,ison2}];
            else
                RNEW = RA*A_NEW{BOX.RHO,ison}; 
            end
        end
    end
    
    if ibox>1
        % fixed skeleton size
        k = size(RNEW,1); 
        % ID given known rank k (fast)
        [T,J] = ID(RNEW,k);
        
        % Record the outgoing skeletons:
        A_NEW{BOX.I_SKUP,ibox} = indskel(J(1:k));
        A_NEW{BOX.T_UP,ibox} = T;
        A_NEW{BOX.J_UP,ibox} = J;
        
        A_NEW{BOX.NSKEL,ibox} = length(indskel);
        A_NEW{BOX.KSKEL,ibox} = k; 
        
        % Record pieces of R and L that need to be merged on the next level
        % rho = R(:,Jup(1:k))
        A_NEW{BOX.RHO,ibox} = RNEW(:,J(1:k));
    
        % Clear auxiliary data
        A_NEW{36,ibox} = [];
        A_NEW{37,ibox} = []; 
    end
end