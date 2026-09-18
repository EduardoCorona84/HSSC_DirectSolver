function T = interpolation_operator(cB,lev,par,X_proxy,X1,X2,type,lr,q)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%     FUNCTION CALL:
%         T = interpolation_operator(cB,lev,par,X_proxy,X1,X2,type,lr,q)
%
%     DESCRIPTION:
%         This function calculates the interpolation operators T that comprise L and R in the
%         telescoping factorization. To speedup computation we use equivalent densities, and if
%         lr=true, we compress K(X_proxy,X2) as low-rank by a randomized ID. We solve for T by
%         least squares:
%
%                     T(X1,X2) = K(X_proxy,X1) \ K(X_proxy,X2)
%
%         In other words, for
%
%                     phi_1 = T(X1,X2) phi_2
%                     K(X_proxy,X1) phi_1 = K(X_proxy,X2) phi_2
%
%         This function corresponds to Algorithm 3 in the paper without the HSS1D compression.
%         It is called within HSS2D_bintree_* during the BUILD TREE stage for boxes with
%         k < par.n_cut.
%
%     INPUT:
%         cB      <2x1 float>         Center of the box.
%         lev     <int>               Level in the tree.
%         par     <struct>            (see HSS_tree_parameters.m) Solve parameters.
%         X_proxy <mx2 float>         Proxy points used in place of X_ext.
%         X1      <kx2 float>         Skeleton points.
%         X2,     <(n-k)x2 float>        Residual points.
%         type    <string>            {'box_2_proxy', 'proxy_2_box'} Direction of interpolation.
%         lr      <bool>              Whether T is stored densely or as low-rank. If true,
%                                         K(X_proxy,X2) is compressed as low-rank using a random
%                                         ID (LRID_rand.m).
%         q       <(depth+1)x1 int>   Array of interpolation operator ranks (used in randomized
%                                         ID to produce a first guess).
%
%     OUTPUT:
%         T[X1,X2] <kx(n-k) float>    Interpolation Operator stored densely or as a low-rank UV'
%                                     pair. If q is the rank, then U is kxq and V is qxm
%




acc = par.acc; 

% Check to see if the inputs are non-empty
if (size(X1,1) > 0 && size(X2,1) > 0)    
    %T is then K(X_proxy,X1)^(-1) * K(X_proxy,X2)
    
    % If lr is 0, we build the operator densely. 
    if lr==0
        if strcmp(type,'box_2_proxy')
            % K[X_ext,X1]
            K11 = Kernel_Eval(X_proxy,X1,par); 
            % K_[X_ext,X2]
            K21 = Kernel_Eval(X_proxy,X2,par); 
            clear X; 
            
            % K[X_ext,X1] x T = K[X_ext,X2]
            T = K11\K21;
        else
            if strcmp(type,'proxy_2_box')
            % K[X1,X_ext]
            K11 = Kernel_Eval(X1,X_proxy,par); 
            % K_[X2,X_ext]
            K21 = Kernel_Eval(X2,X_proxy,par); 
            clear X;     
            
            % K^T[X_ext,X1] x T = K^T[X_ext,X2]
            T = K11.'\K21.'; 
            end
        end
    else
        
        % Else, we compress K21 and use this to speed up the solve (since
        % it is low rank)
        if strcmp(type,'box_2_proxy')
                % K[X_ext,X1]
                K11 = Kernel_Eval(X_proxy,X1,par); 
                % K_[X_ext,X2]
                K21 = Kernel_Eval(X_proxy,X2,par); 
                clear X; 
            
                % Low Rank decomposition using Randomized ID and dense
                % matvecs (HSS1D compression and matvecs is possible, but 
                % more expensive for the sizes considered): 
                mvK = @(x) K21*x; mvKt = @(x) K21.'*x;
                m = size(K21,1); n = size(K21,2); 
                
                % Randomized ID parameters. q is an array storing ranks of
                % interpolation operators in finer levels. 
                
                % If we already have data for ranks on this level, we use
                % that as a first guess. 
                if q(lev+1)>0
                    ID_params.q = q(lev+1)+1; 
                    ID_params.C = 10;
                else
                    % Else, we use ranks on finer levels and a logarithmic
                    % growth model
                    if sum(q>0) < 5
                        % If we don't have enough data, we make a default guess
                        ID_params.q = min(ceil(9*log2(min(m,n))-20),min(m,n));
                        ID_params.C = 5; 
                    else
                        % q(lev+1) ~ C*log2(N/2^lev)+D implies 
                        % q(lev-1) ~ q(lev+1) + 2C ~ q(lev-1) + (q(lev-1)-q(lev-3)).  
                        ID_params.q = min(max(ceil(2*q(lev+3) - q(lev+5) + 2),q(lev+3)),min(m,n)); 
                        ID_params.C = ceil(abs(q(lev+3) - q(lev+5))/2) + 1;
                    end
                end
                
                ID_params.par = acc;  ID_params.form = 'dense'; 
                
                % Lowrank compression using Randomized ID 
                K21_LR = LRID_rand(m,n,mvK,mvKt,'norm',ID_params);
                %K21_LR = LRID(K21,1e-10);
                
                % T = K11\K21
                T.U = K11\K21_LR.U; 
                T.V = K21_LR.V;
        else
                % K[X1,X_ext]
                K11 = Kernel_Eval(X1,X_proxy,par); 
                % K_[X2,X_ext]
                K21 = Kernel_Eval(X2,X_proxy,par); 
                clear X;
                 
                % Low Rank decomposition using Randomized ID and dense
                % matvecs: 
                mvK = @(x) K21.'*x; mvKt = @(x) K21*x;
                m = size(K21',1); n = size(K21',2); 
                
                % Randomized ID parameters. q is an array storing ranks of
                % interpolation operators in finer levels. 
                if q(lev+1)>0
                    ID_params.q = q(lev+1)+1; 
                    ID_params.C = 10;
                else
                    if sum(q>0) < 5
                        ID_params.q = min(ceil(9*log2(min(m,n))-20),min(m,n));
                        ID_params.C = 5; 
                    else
                        ID_params.q = min(ceil(2*q(lev+3) - q(lev+5) + 2),min(m,n)); 
                        ID_params.C = ceil(q(lev+3) - q(lev+5)/2) + 1;
                    end
                end
                
                ID_params.par = acc; ID_params.form = 'dense'; 
                
                K21_LR = LRID_rand(m,n,mvK,mvKt,'norm',ID_params);
                %K21_LR = LRID(K21',1e-10);
                T.U = K11.'\K21_LR.U; 
                T.V = K21_LR.V;
        end
    end

else
    
    % Returns empty array
    T = []; 

end

end   