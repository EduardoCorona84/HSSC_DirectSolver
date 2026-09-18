function A_NEW = HSS1D_recompress_nsym(A,acc)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        A_NEW = HSS1D_recompress_nsym(A,acc)
%
%    DESCRIPTION:
%        This function that recompresses a non-symmetric HSS matrix A, obtaining an optimal
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
            
            L1 = [eye(k1) A{BOX.T_DN,ison1}].'; L1(A{BOX.J_DN,ison1},:) = L1;
            L2 = [eye(k2) A{BOX.T_DN,ison2}].'; L2(A{BOX.J_DN,ison2},:) = L2;
           
            B1 = A{BOX.M_SIB,ison1}; 
            B2 = A{BOX.M_SIB,ison2};
            
            if ibox==1 
                [P1,D1,Q2] = svd(B1,0);
                [P2,D2,Q1] = svd(B2,0);
                
                k1_NEW = sum(diag(D1)>acc); k2_NEW = sum(diag(D2)>acc);
                
                P1 = P1(:,1:k1_NEW); D1 = D1(1:k1_NEW,1:k2_NEW); Q2 = Q2(:,1:k2_NEW); 
                P2 = P2(:,1:k2_NEW); D2 = D2(1:k2_NEW,1:k1_NEW); Q1 = Q1(:,1:k1_NEW);
                
                % Sibling interactions 
                A_NEW{BOX.M_SIB,ison2} = D2; A_NEW{BOX.M_SIB,ison1} = D1;
                % \tilde{T}_i
                A_NEW{36,ison1} = D2; A_NEW{36,ison2} = D1;
                % \tilde{S}_i
                A_NEW{37,ison1} = D1; A_NEW{37,ison2} = D2;
                % \hat{R}_i 
                A_NEW{38,ison1} = Q1'*R1; 
                A_NEW{38,ison2} = Q2'*R2;
                % \hat{L}_i
                A_NEW{39,ison1} = L1*P1; 
                A_NEW{39,ison2} = L2*P2;
            
            else
               % \tilde{T} for parent and R_hat
               Tp = A_NEW{36,ibox};
               Sp = A_NEW{37,ibox};   
               Rp = A_NEW{38,ibox};
               Lp = A_NEW{39,ibox};
               
               % Define \tilde{T} for children 
               Sson1 = [B1   Lp(1:k1,:)*Sp];
               Sson2 = [B2   Lp(k1+1:end,:)*Sp];
               
               Tson2 = [B1 ; Tp*Rp(:,k1+1:end)];
               Tson1 = [B2 ; Tp*Rp(:,1:k1)];
               
               % clear Tp and Sp
               A_NEW{36,ibox} = []; A_NEW{37,ibox} = [];
               
               % SVD(T_i)
                [Y2,DT2,Q2] = svd(Tson2,0); k2_NEW = sum(diag(DT2)>acc); 
                [Y1,DT1,Q1] = svd(Tson1,0); k1_NEW = sum(diag(DT1)>acc); 
                
               % SVD(S_i)
                [P1,DS1,~] = svd(Sson1,0); k1_NEW = max(sum(diag(DS1)>acc),k1_NEW); 
                [P2,DS2,~] = svd(Sson2,0); k2_NEW = max(sum(diag(DS2)>acc),k2_NEW);
            
                % Re-size the factors of the svd
                Y1 = Y1(:,1:k1_NEW); DT1 = DT1(1:k1_NEW,1:k1_NEW); Q1 = Q1(:,1:k1_NEW);
                Y2 = Y2(:,1:k2_NEW); DT2 = DT2(1:k2_NEW,1:k2_NEW); Q2 = Q2(:,1:k2_NEW);
                P1 = P1(:,1:k1_NEW); DS1 = DS1(1:k1_NEW,1:k1_NEW); 
                P2 = P2(:,1:k2_NEW); DS2 = DS2(1:k2_NEW,1:k2_NEW); 
            
                % Sibling Interactions
                A_NEW{BOX.M_SIB,ison1} = P1'*Y2(1:k1,:)*DT2; 
                A_NEW{BOX.M_SIB,ison2} = P2'*Y1(1:k2,:)*DT1; 
                % \tilde{T} (diagonal matrix)
                A_NEW{36,ison1} = DT1; A_NEW{36,ison2} = DT2;
                % \tilde{S} (diagonal matrix)
                A_NEW{37,ison1} = DS1; A_NEW{37,ison2} = DS2;
                % \hat{R}_i 
                A_NEW{38,ison1} = Q1'*R1; 
                A_NEW{38,ison2} = Q2'*R2;
                % \hat{L}_i
                A_NEW{39,ison1} = L1*P1; 
                A_NEW{39,ison2} = L2*P2;
                
                % \tilde{R}_p ,\tilde{L}_p (new "interpolation" operators)
                A_NEW{38,ibox} = [Rp(:,1:k1)*Q1 Rp(:,k1+1:end)*Q2]; 
                A_NEW{39,ibox} = [P1'*Lp(1:k1,:) ; P2'*Lp(k1+1:end,:)];
            end
            
        else
            ison = max(ison1,ison2);
            k = A{BOX.KSKEL,ison}; 
            R = [eye(k) A{BOX.T_UP,ison}]; R(:,A{BOX.J_UP,ison}) = R;
            L = [eye(k) A{BOX.T_DN,ison}].'; L(A{BOX.J_DN,ison},:) = L;
            
            Tp = A_NEW{36,ibox};
            Sp = A_NEW{37,ibox};
            Rp = A_NEW{38,ibox};
            Lp = A_NEW{39,ibox};
            
            % Since there are no sibling interactions, this box interacts
            % with the neutered column/row only through its parent. 
            Sson = Lp*Sp;
            Tson = Tp*Rp;
            
            % clear Tp and Sp
            A_NEW{36,ibox} = []; A_NEW{37,ibox} = [];
            
            % SVD(T_i)
            [~,DT,Q] = svd(Tson,0); k_NEW = sum(diag(DT)>acc);
            [P,DS,~] = svd(Sson,0); k_NEW = max(sum(diag(DS)>acc),k_NEW);
            % Re-size the factors of the svd
            DT = DT(1:k_NEW,1:k_NEW); Q = Q(:,1:k_NEW);
            DS = DS(1:k_NEW,1:k_NEW); P = P(:,1:k_NEW);
            
            % \tilde{T} 
            A_NEW{36,ison} = DT;
            % \tilde{S} 
            A_NEW{37,ison} = DS;
            % \hat{R}_i 
            A_NEW{38,ison} = Q'*R;
            % \hat{L}_i 
            A_NEW{39,ison} = L*P;
            % \tilde{R}_p , \tilde{L}_p 
            A_NEW{38,ibox} = Rp*Q; 
            A_NEW{39,ibox} = P'*Lp;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (II) Re-formatting to ID HSS: 

for ibox = nboxes:(-1):1
    if ( (A{BOX.C1,ibox}<=0) && (A{BOX.C2,ibox}<=0) ) % ibox has no sons.
        %indices
        indskel_out = A{BOX.END1,ibox} - 1 + (1:A{BOX.END2,ibox});
        indskel_in  = indskel_out; 
        % Self-interactions are the same
        A_NEW{BOX.M_SELF,ibox} = A{BOX.M_SELF,ibox};
        % We need to re-format the new interpolation operator
        RNEW = A_NEW{38,ibox};
        LNEW = A_NEW{39,ibox};
    else
        ison1       = A_NEW{BOX.C1,ibox};
        ison2       = A_NEW{BOX.C2,ibox};
        if (ison1>0 && ison2>0)
            indskel_out = [A_NEW{BOX.I_SKUP,ison1},A_NEW{BOX.I_SKUP,ison2}];
            indskel_in  = [A_NEW{BOX.I_SKDN,ison1},A_NEW{BOX.I_SKDN,ison2}];
            
             % Sibling Interaction Matrices
            A_NEW{BOX.M_SIB,ison1} = A_NEW{BOX.LAMBDA,ison1}*A_NEW{BOX.M_SIB,ison1}*A_NEW{BOX.RHO,ison2};
            A_NEW{BOX.M_SIB,ison2} = A_NEW{BOX.LAMBDA,ison2}*A_NEW{BOX.M_SIB,ison2}*A_NEW{BOX.RHO,ison1}; 
        else
            ison = max([ison1,ison2]); 
            indskel_out = A_NEW{BOX.I_SKUP,ison};
            indskel_in  = A_NEW{BOX.I_SKDN,ison};
        end
        
        if ibox>1
        
            % Interpolation Matrices 
            RA = A_NEW{38,ibox};
            LA = A_NEW{39,ibox};
            
            % We include factors coming from recompression of leaves. 
            if (ison1>0 && ison2>0)
                kA1 = size(A_NEW{BOX.RHO,ison1},1); 
                RNEW = [RA(:,1:kA1)*A_NEW{BOX.RHO,ison1} RA(:,kA1+1:end)*A_NEW{BOX.RHO,ison2}];
                LNEW = [A_NEW{BOX.LAMBDA,ison1}*LA(1:kA1,:) ; A_NEW{BOX.LAMBDA,ison2}*LA(kA1+1:end,:)];
            else
                RNEW = RA*A_NEW{BOX.RHO,ison}; 
                LNEW = A_NEW{BOX.LAMBDA,ison}*LA;
            end
        end
    end
    
    if ibox>1
        % fixed skeleton size (from svds on previous step)
        k = min(size(RNEW,1),size(RNEW,2)); 
        % ID given known rank k (fast)
        [Tout,Jout] = ID(RNEW,k);
        [Tin,Jin] = ID(LNEW.',k);
        
        % Record the outgoing skeletons:
        A_NEW{BOX.I_SKUP,ibox} = indskel_out(Jout(1:k));
        A_NEW{BOX.T_UP,ibox} = Tout;
        A_NEW{BOX.J_UP,ibox} = Jout;
        
        % Record the outgoing skeletons:
        A_NEW{BOX.I_SKDN,ibox} = indskel_in(Jin(1:k));
        A_NEW{BOX.T_DN,ibox} = Tin;
        A_NEW{BOX.J_DN,ibox} = Jin;
        
        A_NEW{BOX.NSKEL,ibox} = length(indskel_out);
        A_NEW{BOX.KSKEL,ibox} = k; 
        
        % Record pieces of R and L that need to be merged on the next level
        % rho = R(:,Jup(1:k))
        A_NEW{BOX.RHO,ibox} = RNEW(:,Jout(1:k));
        A_NEW{BOX.LAMBDA,ibox} = LNEW(Jin(1:k),:);
    
        % Clear auxiliary data
        A_NEW{36,ibox} = [];
        A_NEW{37,ibox} = []; 
    end
end