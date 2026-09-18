function x = HSS1D_OLS(KHSS,b,acc,tau,mu,maxit,type)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        x = HSS1D_OLS(KHSS,b,acc,tau,mu,maxit,type)
%
%    DESCRIPTION:
%        This function solves the system Ax=b for x using the HSS1D overdetermined least squares
%        algorithm based on Barlow & Vemulapati's approach (Ken & Greengard 12'). Alternatively,
%        one can compute using the corrected semi-normal equations R'RX = A'b where A = QR.
%
%        (1) The HSS1D structure for A is turned into an augmented sparse system using equivalent
%            densities as auxiliary variables.
%        (2) The least squares problem for A turns into a linearly constrained least squares of the
%            form min{|Ex-b|} s.t. Cx=0
%        (3) Barlow & Vemulapatti's approach is to turn this into an iteration performing iterative
%            refinement on a non-constrained, weighted least squares of the form
%            min{|Ex-b|+tau|Cx|}
%
%    INPUT:
%        KHSS is the {50 x nboxes} cell array that encodes matrix A in HSS form.
%        b is the right hand side (N x 1) vector.
%        acc (double) is the desired accuracy.
%        tau is a positive weight on linear constraints. The optimal choice for tau for precision
%            eps (1e-16 for double) is tau=(mu*cond(A))^(-1/3)
%        mu is a small, positive number for Tikhonov regularization. If mu=0, no regularization is
%            performed.
%        maxit max number of iterations
%        type is a string with options 'qr' to use QR decomposition for the iteration, and 'sn' to
%            use semi-normal equations.
%
%    OUTPUT:
%        x is the least squares solution vector of Ax = b.
%


global BOX

% Create augmented sparse matrix and multiply auxiliary sub-block by tau
display(log10(mu))
Ksp = HSS1D_sparse_matrix(KHSS,tau); 
[Msp,Nsp] = size(Ksp); 


% Dimensions of original matrix
N = KHSS{BOX.END2,1};
% Check if matrix is rectangular. For now we assume M>N. 
rect = ~isempty(KHSS{17,1}); 
if rect
    M = KHSS{17,1}; 
else
    M = N; 
end

if mu>0
   [I,J,K] = find(Ksp); 
   ind_var = I<M+1;
   I_var = I(ind_var);  J_var = J(ind_var); K_var = K(ind_var); 
   I_reg = (M+1:M+N).'; J_reg = (1:N).';    K_reg = mu*ones(N,1); 
   I_ctr = I(~ind_var)+N; J_ctr = J(~ind_var); K_ctr = K(~ind_var); 
   
   Ksp = sparse([I_var;I_reg;I_ctr],[J_var;J_reg;J_ctr],[K_var;K_reg;K_ctr],Msp+N,Nsp); 
   Msp = Msp+N; 
   
   ind_var = 1:M; 
   ind_ctr = M+N+1:Msp; 
   ind_reg = M+1:M+N; 
else
    ind_var = 1:M; 
    ind_ctr = M+1:Msp; 
end

B = zeros(Msp,size(b,2));
B(ind_var,:) = b; 

if strcmp(type,'qr')
    %E = sparse(I_var,J_var,K_var,M,Nsp);
    %C = sparse(I_ctr-(M+N),J_ctr,K_ctr,Msp-(M+N),Nsp); 
    
    %fprintf('\n Condition Number \n')
    %display(cond(full(Ksp)));          
    
    %sigma = gsvd(full(C),full(E));  
    %fprintf('\n Generalized Singular Values \n')
    %display(max(sigma(sigma<Inf)));     
    
    % Economy size QR. Returns C = Q'*B and sparse, triangular R. 
    [C,R] = qr(Ksp,B,0); 
    X = R\C; 

    %Residual
    Res = B - Ksp*X;
    lambda = tau*Res(ind_ctr); 

    Err = norm(Res(ind_ctr))/tau;  
    iter = 1; 
    display(iter)
    display(norm(Res(ind_var)))
    display(norm(Res(ind_ctr)))

    while Err>acc && iter<maxit
        Bk = Res; 
        Bk(ind_ctr) = Bk(ind_ctr) + (1/tau)*lambda; 
    
        % Here we should really be re-using the sparse QR. Consider solving
        % semi-normal equations instead. 
        
        if mu>0
           [I,J,K] = find(Ksp); 
           K(I>M & I<M+N+1) = 0.5*K(I>M & I<M+N+1);
           Ksp = sparse(I,J,K,Msp,Nsp); 
        end
        
        [C,R] = qr(Ksp,Bk,0); 
        dX = R\C; 
    
        X = X + dX;
        Res = Res - Ksp*dX;
        lambda = tau*Res(ind_ctr); 
        Err = norm(Res(ind_ctr))/tau;  
        iter = iter+1; 
        
        display(iter)
        display(norm(Res(ind_var)))
        display(norm(Res(ind_ctr)))
    end

elseif strcmp(type,'sn')
    % Solve w/ semi-normal equations and iterative refinement
    % Economy size QR. Returns sparse, triangular R. 
    R = qr(Ksp,0);
    
    %X = (R.'*R + mu*eye(size(R)))\(Ksp.'*B); 
    X = R\(R.'\(Ksp.'*B));
    
    %Residual
    Res = B - Ksp*X;

    Err = norm(Res(ind_ctr))/tau; 
    display(norm(Res(ind_var)))
    display(norm(Res(ind_ctr)))
    iter = 1; 
    display(iter)
    
    while Err>acc && iter<maxit
        
        %dX = (R.'*R + mu*eye(size(R)))\(Ksp.'*Res); 
        dX = R\(R.'\(Ksp.'*Res));
    
        X = X + dX;
        Res = Res - Ksp*dX;
        Err = norm(Res(ind_ctr))/tau;
        display(norm(Res(ind_var)))
        display(norm(Res(ind_ctr)))
        iter = iter+1; 
        display(iter)
    end
end

x = X(1:N,:); 
%display(Err)
%display(iter)
