function KHSS_sp = HSS1D_sparse_matrix(KHSS,tau)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        KHSS_sp = HSS1D_sparse_matrix(KHSS,tau)
%
%    DESCRIPTION:
%        Given an NxN matrix in HSS1D form, this code produces the corresponding augmented sparse
%        system of size S ~ NlogN. This can then be used to compute fast arithmetic of HSS matrices
%        or perform sparse linear algebra such as sparse QR, Cholesky, Least Squares, etc.
%
%    INPUT:
%        KHSS is the HSS structure of a matrix.
%        tau (double) is an option for the OLS code. The submatrix corresponding to rows of
%        auxiliary variables (upward and downward 'equivalent densities') is multiplied by tau.
%
%    OUTPUT:
%        KHSS_sp (sparse, block-banded matrix) =
%                [D^d L^d  0 ...                   |
%                |R^d  0  -I ...                   |
%                |    -I  D^{d-1} ...              |
%                |    ...              -I          |
%                |                 -I   D^1 L^1    |
%                |                      R^1  0  -I |
%                |                          -I  D^0]
%


global BOX

nboxes = size(KHSS,2); 
depth = KHSS{BOX.LEVEL,nboxes}; %Tree depth

% Check if matrix is rectangular or not
rect = ~isempty(KHSS{17,1}); 

% Count of box sizes per level
N = zeros(depth,1); Sk = N; 
if rect
    M = N; Sk_in = Sk; 
end

for i=2:nboxes                                        

% Auxiliary indices
if N(KHSS{BOX.LEVEL,i}) == 0
    KHSS{10,i} = 1; 
    KHSS{11,i} = 1; 
    if rect
       KHSS{12,i} = 1; 
       KHSS{13,i} = 1; 
    end
else
    KHSS{10,i} = KHSS{10,i-1} + KHSS{BOX.NSKEL,i-1}; 
    KHSS{11,i} = KHSS{11,i-1} + KHSS{BOX.KSKEL,i-1};
    if rect
       KHSS{12,i} = KHSS{12,i-1} + KHSS{18,i-1}; 
       KHSS{13,i} = KHSS{13,i-1} + KHSS{19,i-1};
    end
end

N(KHSS{BOX.LEVEL,i})  = N(KHSS{BOX.LEVEL,i}) + KHSS{BOX.NSKEL,i}; 
Sk(KHSS{BOX.LEVEL,i}) = Sk(KHSS{BOX.LEVEL,i}) + KHSS{BOX.KSKEL,i};
if rect
    M(KHSS{BOX.LEVEL,i})     = M(KHSS{BOX.LEVEL,i})     + KHSS{18,i}; 
    Sk_in(KHSS{BOX.LEVEL,i}) = Sk_in(KHSS{BOX.LEVEL,i}) + KHSS{19,i};
end

end

KHSS{10,1} = 1; 
KHSS{11,1} = 1; 
KHSS{BOX.NSKEL,1} = KHSS{BOX.KSKEL,2}+KHSS{BOX.KSKEL,3}; 

if rect
    KHSS{12,1} = 1; 
    KHSS{13,1} = 1; 
    KHSS{18,1} = KHSS{19,2}+KHSS{19,3}; 
end

%KHSS(1:19,:)

% Sparse matrix size
S = 0; 
if rect
   Sm = 0;  
end

% KHSS_sp will be = sparse(I,J,K,M,M). In other words, the sparse structure
% is represented by vectors I,J and K such that K(I,J) = M and only nonzero entries 
% are stored. 
I = []; 
J = []; 
K = []; 

for ibox = nboxes:(-1):1
    
    if ( (KHSS{BOX.C1,ibox}<=0) && (KHSS{BOX.C2,ibox}<=0) ) % ibox has no sons.
        
        if rect
            indN_out = KHSS{10,ibox} - 1 + (1:KHSS{BOX.NSKEL,ibox});
            indN_in  = KHSS{12,ibox} - 1 + (1:KHSS{18,ibox});
            indK_out = M(KHSS{BOX.LEVEL,ibox}) + KHSS{11,ibox} - 1 + (1:KHSS{BOX.KSKEL,ibox});
            indK_in  = N(KHSS{BOX.LEVEL,ibox}) + KHSS{13,ibox} - 1 + (1:KHSS{19,ibox});
        
            % Create meshgrids of indices for D,L and R
            [ID,JD] = meshgrid(indN_in,indN_out); 
            [IL,JL] = meshgrid(indN_in,indK_in);
            [IR,JR] = meshgrid(indK_out,indN_out);  
        
            % Compute interpolation matrix
            n = KHSS{BOX.NSKEL,ibox};  k_out = KHSS{BOX.KSKEL,ibox}; 
            m = KHSS{18,ibox}; k_in  = KHSS{19,ibox}; 
            Jout = KHSS{BOX.J_UP,ibox}; R = [eye(k_out) KHSS{BOX.T_UP,ibox}];     
            R(:,Jout) = R; 
            Jin = KHSS{BOX.J_DN,ibox}; L = [eye(k_in) KHSS{BOX.T_DN,ibox}].';     
            L(Jin,:) = L; 
       
            % Augment index vectors (I,J) and entry vector K. This adds 
            % ... D^d L^d  0 ...                   
            % ... R^d  0  -I ...                   
            %     -I  
            % to matrix KHSS_sp. 
            
            I = [I ; ID(:) ; IL(:) ; IR(:) ; indK_out' ; indK_in' + (M(KHSS{BOX.LEVEL,ibox}) - N(KHSS{BOX.LEVEL,ibox})) + Sk(KHSS{BOX.LEVEL,ibox})]; 
            J = [J ; JD(:) ; JL(:) ; JR(:) ; indK_out' + (N(KHSS{BOX.LEVEL,ibox}) - M(KHSS{BOX.LEVEL,ibox})) + Sk_in(KHSS{BOX.LEVEL,ibox}) ; indK_in'];
            K = [K ; reshape(KHSS{BOX.M_SELF,ibox}.',m*n,1) ; ...
                 reshape(L.',m*k_in,1) ; tau*reshape(R.',n*k_out,1) ; ...
                 -tau*ones(k_out,1) ; -tau*ones(k_in,1)]; 
        else
            indN = KHSS{10,ibox} - 1 + (1:KHSS{BOX.NSKEL,ibox});
            indK = N(KHSS{BOX.LEVEL,ibox}) + KHSS{11,ibox} - 1 + (1:KHSS{BOX.KSKEL,ibox});
        
            % Create meshgrids of indices for D,L and R
            [ID,JD] = meshgrid(indN,indN); 
            [IL,JL] = meshgrid(indN,indK);
            [IR,JR] = meshgrid(indK,indN);  
        
            % Compute interpolation matrix
            n = KHSS{BOX.NSKEL,ibox}; k = KHSS{BOX.KSKEL,ibox}; 
            Jout = KHSS{BOX.J_UP,ibox}; R = [eye(k) KHSS{BOX.T_UP,ibox}];     
            R(:,Jout) = R; 
            if isempty(KHSS{BOX.T_DN,ibox})
                L = R.';  
            else
                Jin = KHSS{BOX.J_DN,ibox}; L = [eye(k) KHSS{BOX.T_DN,ibox}].';     
                L(Jin,:) = L;  
            end
       
            % Augment index vectors (I,J) and entry vector K. This adds 
            % ... D^d L^d  0 ...                   
            % ... R^d  0  -I ...                   
            %     -I  
            % to matrix KHSS_sp. 
            I = [I ; ID(:) ; IL(:) ; IR(:) ; indK' ; indK' + Sk(KHSS{BOX.LEVEL,ibox})]; 
            J = [J ; JD(:) ; JL(:) ; JR(:) ; indK' + Sk(KHSS{BOX.LEVEL,ibox}) ; indK'];
            K = [K ; reshape(KHSS{BOX.M_SELF,ibox}.',n^2,1) ; ...
                 reshape(L.',n*k,1) ; tau*reshape(R.',n*k,1) ; ...
                 -tau*ones(k,1) ; -tau*ones(k,1)]; 
        end
             
        if KHSS{10,ibox} == 1
           % Increase matrix length S
          
           if rect
               Sm = Sm + M(KHSS{BOX.LEVEL,ibox}) + Sk(KHSS{BOX.LEVEL,ibox});  
               S  = S +  N(KHSS{BOX.LEVEL,ibox}) + Sk_in(KHSS{BOX.LEVEL,ibox});  
           else
                S = S + N(KHSS{BOX.LEVEL,ibox}) + Sk(KHSS{BOX.LEVEL,ibox});  
           end
        end
    
    elseif ( (KHSS{BOX.C1,ibox}>0) && (KHSS{BOX.C2,ibox}>0) )
        
        ison1       = KHSS{BOX.C1,ibox}; k1 = KHSS{BOX.KSKEL,ison1}; 
        ison2       = KHSS{BOX.C2,ibox}; k2 = KHSS{BOX.KSKEL,ison2}; 
        
        if rect
            
            k1_in = KHSS{19,ison1}; k2_in = KHSS{19,ison2}; 
            D = [zeros(k1_in,k1) tau*KHSS{BOX.M_SIB,ison1} ; tau*KHSS{BOX.M_SIB,ison2} zeros(k2_in,k2)]; 
        
        % If ibox>1, we add matrices D,L and R. If ibox==1, we need only
        % add D^0. 
        if ibox>1
            indN_out = S  + KHSS{10,ibox} - 1 + (1:KHSS{BOX.NSKEL,ibox});
            indN_in  = Sm + KHSS{12,ibox} - 1 + (1:KHSS{18,ibox});
            indK_out = Sm + M(KHSS{BOX.LEVEL,ibox}) + KHSS{11,ibox} - 1 + (1:KHSS{BOX.KSKEL,ibox});
            indK_in  = S  + N(KHSS{BOX.LEVEL,ibox}) + KHSS{13,ibox} - 1 + (1:KHSS{19,ibox});
        
            % Create meshgrids of indices for D,L and R
            [ID,JD] = meshgrid(indN_in,indN_out); 
            [IL,JL] = meshgrid(indN_in,indK_in);
            [IR,JR] = meshgrid(indK_out,indN_out);  
        
            % Compute interpolation matrix
            n = KHSS{BOX.NSKEL,ibox};  k_out = KHSS{BOX.KSKEL,ibox}; 
            m = KHSS{18,ibox}; k_in  = KHSS{19,ibox}; 
            Jout = KHSS{BOX.J_UP,ibox}; R = [eye(k_out) KHSS{BOX.T_UP,ibox}];     
            R(:,Jout) = R; 
            Jin = KHSS{BOX.J_DN,ibox}; L = [eye(k_in) KHSS{BOX.T_DN,ibox}].';     
            L(Jin,:) = L; 
       
            % Augment index vectors (I,J) and entry vector K. This adds 
            % ... D^d L^d  0 ...                   
            % ... R^d  0  -I ...                   
            %     -I  
            % to matrix KHSS_sp. 
            I = [I ; ID(:) ; IL(:) ; IR(:) ; ...
                indK_out' ; indK_in' + (Sm + M(KHSS{BOX.LEVEL,ibox}) - S - N(KHSS{BOX.LEVEL,ibox})) + Sk(KHSS{BOX.LEVEL,ibox})]; 
            J = [J ; JD(:) ; JL(:) ; JR(:) ; ...
                indK_out' + (S + N(KHSS{BOX.LEVEL,ibox}) - Sm - M(KHSS{BOX.LEVEL,ibox})) + Sk_in(KHSS{BOX.LEVEL,ibox}) ; indK_in'];
            K = [K ; reshape(D.',m*n,1) ; ...
                 tau*reshape(L.',m*k_in,1) ; tau*reshape(R.',n*k_out,1) ; ...
                 -tau*ones(k_out,1) ; -tau*ones(k_in,1)]; 
        else
            indN_out = S  + KHSS{10,ibox} - 1 + (1:KHSS{BOX.NSKEL,ibox});
            indN_in  = Sm + KHSS{12,ibox} - 1 + (1:KHSS{18,ibox});           
            [ID,JD] = meshgrid(indN_in,indN_out); 
            n = KHSS{BOX.NSKEL,ibox}; m = KHSS{18,ibox}; 
            
            % Augment index vectors (I,J) and entry vector K.
            I = [I ; ID(:)]; 
            J = [J ; JD(:)];            
            K = [K ; reshape(D.',m*n,1) ]; 
        end
        
        else
        D = [zeros(k1,k1) tau*KHSS{BOX.M_SIB,ison1} ; tau*KHSS{BOX.M_SIB,ison2} zeros(k2,k2)]; 
        
        % If ibox>1, we add matrices D,L and R. If ibox==1, we need only
        % add D^0. 
        if ibox>1
            indN = S + KHSS{10,ibox} - 1 + (1:KHSS{BOX.NSKEL,ibox});
            indK = S + N(KHSS{BOX.LEVEL,ibox}) + KHSS{11,ibox} - 1 + (1:KHSS{BOX.KSKEL,ibox});
            
            % Create meshgrids of indices for D,L and R
            [ID,JD] = meshgrid(indN,indN); 
            [IL,JL] = meshgrid(indN,indK);
            [IR,JR] = meshgrid(indK,indN);  
            
            % Compute interpolation matrix
            n = KHSS{BOX.NSKEL,ibox}; k = KHSS{BOX.KSKEL,ibox}; 
            Jout = KHSS{BOX.J_UP,ibox}; 
            R = [eye(k) KHSS{BOX.T_UP,ibox}];     
            R(:,Jout) = R; 
            if isempty(KHSS{BOX.T_DN,ibox})
                L = R.';  
            else
                Jin = KHSS{BOX.J_DN,ibox}; L = [eye(k) KHSS{BOX.T_DN,ibox}].';     
                L(Jin,:) = L;  
            end
        
            % Augment index vectors (I,J) and entry vector K.
            I = [I ; ID(:) ; IL(:) ; IR(:) ; indK' ; indK' + Sk(KHSS{BOX.LEVEL,ibox})]; 
            J = [J ; JD(:) ; JL(:) ; JR(:) ; indK' + Sk(KHSS{BOX.LEVEL,ibox}) ; indK'];
            K = [K ; reshape(D.',n^2,1) ; ...
                 tau*reshape(L.',n*k,1) ; tau*reshape(R.',n*k,1) ; ...
                 -tau*ones(k,1) ; -tau*ones(k,1)]; 
        else
            indN = S + KHSS{10,ibox} - 1 + (1:KHSS{BOX.NSKEL,ibox});            
            [ID,JD] = meshgrid(indN,indN); 
            n = KHSS{BOX.NSKEL,ibox}; 
            
            % Augment index vectors (I,J) and entry vector K.
            I = [I ; ID(:)]; 
            J = [J ; JD(:)];            
            K = [K ; reshape(D.',n^2,1) ]; 
        end
        end
             
        if KHSS{10,ibox} == 1
            % Increase matrix length S
           if rect
           if ibox>1 
               Sm = Sm + M(KHSS{BOX.LEVEL,ibox}) + Sk(KHSS{BOX.LEVEL,ibox});  
               S  = S +  N(KHSS{BOX.LEVEL,ibox}) + Sk_in(KHSS{BOX.LEVEL,ibox});  
           else
               Sm = Sm + m; 
               S = S + n;  
           end
           else
           if ibox>1 
              S = S + N(KHSS{BOX.LEVEL,ibox}) + Sk(KHSS{BOX.LEVEL,ibox});  
           else
              S = S + n;  
           end
           end
        end        
    end
    
  
end

% Create sparse data structure KHSS_sp.
if rect
    KHSS_sp = sparse(I,J,K,Sm,S); 
else
    KHSS_sp = sparse(I,J,K,S,S); 
end

