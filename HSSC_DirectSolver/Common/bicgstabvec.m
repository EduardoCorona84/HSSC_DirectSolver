function [x,flag,relres,iter,resvec] = bicgstabvec(A,b,tol,maxit,M1,M2,x0,varargin)
%BICGSTAB   BiConjugate Gradients Stabilized Method.
%   X = BICGSTAB(A,B) attempts to solve the system of linear equations
%   A*X=B for X. The N-by-N coefficient matrix A must be square and the
%   right hand side column vector B must have length N.
%
%   X = BICGSTAB(AFUN,B) accepts a function handle AFUN instead of the
%   matrix A. AFUN(X) accepts a vector input X and returns the
%   matrix-vector product A*X. In all of the following syntaxes, you can
%   replace A by AFUN.
%
%   X = BICGSTAB(A,B,TOL) specifies the tolerance of the method. If TOL is
%   [] then BICGSTAB uses the default, 1e-6.
%
%   X = BICGSTAB(A,B,TOL,MAXIT) specifies the maximum number of iterations.
%   If MAXIT is [] then BICGSTAB uses the default, min(N,20).
%
%   X = BICGSTAB(A,B,TOL,MAXIT,M) and X = BICGSTAB(A,B,TOL,MAXIT,M1,M2) use
%   preconditioner M or M=M1*M2 and effectively solve the system
%   A*inv(M)*X = B for X. If M is [] then a preconditioner is not
%   applied. M may be a function handle returning M\X.
%
%   X = BICGSTAB(A,B,TOL,MAXIT,M1,M2,X0) specifies the initial guess.  If
%   X0 is [] then BICGSTAB uses the default, an all zero vector.
%
%   [X,FLAG] = BICGSTAB(A,B,...) also returns a convergence FLAG:
%    0 BICGSTAB converged to the desired tolerance TOL within MAXIT iterations.
%    1 BICGSTAB iterated MAXIT times but did not converge.
%    2 preconditioner M was ill-conditioned.
%    3 BICGSTAB stagnated (two consecutive iterates were the same).
%    4 one of the scalar quantities calculated during BICGSTAB became
%      too small or too large to continue computing.
%
%   [X,FLAG,RELRES] = BICGSTAB(A,B,...) also returns the relative residual
%   NORM(B-A*X)/NORM(B). If FLAG is 0, then RELRES <= TOL.
%
%   [X,FLAG,RELRES,ITER] = BICGSTAB(A,B,...) also returns the iteration
%   number at which X was computed: 0 <= ITER <= MAXIT. ITER may be an
%   integer + 0.5, indicating convergence half way through an iteration.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = BICGSTAB(A,B,...) also returns a vector
%   of the residual norms at each half iteration, including NORM(B-A*X0).
%
%   Example:
%      n = 21; A = gallery('wilk',n);  b = sum(A,2);
%      tol = 1e-12;  maxit = 15; M = diag([10:-1:1 1 1:10]);
%      x = bicgstab(A,b,tol,maxit,M);
%   Or, use this matrix-vector product function
%      %-----------------------------------------------------------------%
%      function y = afun(x,n)
%      y = [0; x(1:n-1)] + [((n-1)/2:-1:0)'; (1:(n-1)/2)'].*x+[x(2:n); 0];
%      %-----------------------------------------------------------------%
%   and this preconditioner backsolve function
%      %------------------------------------------%
%      function y = mfun(r,n)
%      y = r ./ [((n-1)/2:-1:1)'; 1; (1:(n-1)/2)'];
%      %------------------------------------------%
%   as inputs to BICGSTAB:
%      x1 = bicgstab(@(x)afun(x,n),b,tol,maxit,@(x)mfun(x,n));
%
%   Class support for inputs A,B,M1,M2,X0 and the output of AFUN:
%      float: double
%
%   See also BICG, BICGSTABL, CGS, GMRES, LSQR, MINRES, PCG, QMR, SYMMLQ,
%   TFQMR, ILU, FUNCTION_HANDLE.

%   Copyright 1984-2010 The MathWorks, Inc.
%   $Revision: 1.19.4.7 $ $Date: 2010/11/17 11:29:41 $

% Check for an acceptable number of input arguments
if nargin < 2
    error(message('MATLAB:bicgstab:NotEnoughInputs'));
end

% Determine whether A is a matrix or a function.
[atype,afun,afcnstr] = iterchk(A);
if strcmp(atype,'matrix')
    % Check matrix and right hand side vector inputs have appropriate sizes
    [m,n] = size(A);
    if (m ~= n)
        error(message('MATLAB:bicgstab:NonSquareMatrix'));
    end
    if ~isequal(size(b),[m,1])
        error(message('MATLAB:bicgstab:RSHsizeMatchCoeffMatrix', m));
    end
else
    m = size(b,1);
    n = m;
    %if ~iscolumn(b)
    %    error(message('MATLAB:bicgstab:RSHnotColumn'));
    %end
end

% Assign default values to unspecified parameters
if nargin < 3 || isempty(tol)
    tol = 1e-6;
end
warned = 0;
if tol < eps
    warning('MATLAB:bicgstab:tooSmallTolerance', ...
        'Input tol is smaller than eps and may not be achieved by BICGSTAB\n         Try to use a bigger tolerance');
    warned = 1;
    tol = eps;
elseif tol >= 1
    warning(message('MATLAB:bicgstab:tooBigTolerance'));
    warned = 1;
    tol = 1-eps;
end
if nargin < 4 || isempty(maxit)
    maxit = min(n,20);
end

% Check for all zero right hand side vector => all zero solution
n2b = norm(b);                     % Norm of rhs vector, b
ns = size(b,2); 
if (n2b == 0)                      % if    rhs vector is all zeros
    x = zeros(n,ns);                % then  solution is all zeros
    flag = 0;                      % a valid solution has been obtained
    relres = 0;                    % the relative residual is actually 0/0
    iter = 0;                      % no iterations need be performed
    resvec = 0;                    % resvec(1) = norm(b-A*x) = norm(0)
    if (nargout < 2)
        itermsg('bicgstab',tol,maxit,0,flag,iter,NaN);
    end
    return
end

if ((nargin >= 5) && ~isempty(M1))
    existM1 = 1;
    [m1type,m1fun,m1fcnstr] = iterchk(M1);
    if strcmp(m1type,'matrix')
        if ~isequal(size(M1),[m,m])
            error(message('MATLAB:bicgstab:WrongPrecondSize', m));
        end
    end
else
    existM1 = 0;
    m1type = 'matrix';
end

if ((nargin >= 6) && ~isempty(M2))
    existM2 = 1;
    [m2type,m2fun,m2fcnstr] = iterchk(M2);
    if strcmp(m2type,'matrix')
        if ~isequal(size(M2),[m,m])
            error(message('MATLAB:bicgstab:WrongPrecondSize', m));
        end
    end
else
    existM2 = 0;
    m2type = 'matrix';
end

if ((nargin >= 7) && ~isempty(x0))
    if ~isequal(size(x0),[n,1])
        error(message('MATLAB:bicgstab:WrongInitGuessSize', n));
    else
        x = x0;
    end
else
    x = zeros(n,ns);
end

if ((nargin > 7) && strcmp(atype,'matrix') && ...
        strcmp(m1type,'matrix') && strcmp(m2type,'matrix'))
    error(message('MATLAB:bicgstab:TooManyInputs'));
end

% Set up for the method
flag = 1;
xmin = x;                          % Iterate which has minimal residual so far
imin = 0;                          % Iteration at which xmin was computed
tolb = tol * n2b;                  % Relative tolerance
r = b - iterapp('mtimes',afun,atype,afcnstr,x,varargin{:});
normr = norm(r);                   % Norm of residual
normr_act = normr;

if (normr <= tolb)                 % Initial guess is a good enough solution
    flag = 0;
    relres = normr / n2b;
    iter = 0;
    resvec = normr;
    if (nargout < 2)
        itermsg('bicgstab',tol,maxit,0,flag,iter,relres);
    end
    return
end

rt = r;                            % Shadow residual
resvec = zeros(2*maxit+1,1);       % Preallocate vector for norm of residuals
resvec(1) = normr;                 % resvec(1) = norm(b-A*x0)
normrmin = normr;                  % Norm of residual from xmin
rho = ones(1,ns);
omega = 1;
stag = 0;                          % stagnation of the method
alpha = [];                        % overshadow any functions named alpha
moresteps = 0;
maxmsteps = min([floor(n/50),10,n-maxit]);
maxstagsteps = 3;
% loop over maxit iterations (unless convergence or failure)

for ii = 1 : maxit
    rho1 = rho; 
    for j=1:ns
    rho(j) = rt(:,j)' * r(:,j);
        if (rho(j) == 0.0) || isinf(rho(j))
            flag = 4;
            resvec = resvec(1:2*ii-1);
            break
        end
    end
    
    if ii == 1
        p = r;
    else
        beta = (rho./rho1).*(alpha./omega);
        if (norm(beta) == 0) | ~isfinite(beta)
            flag = 4;
            break
        end
        
        p = r + (p - v * diag(omega)) * diag(beta);
    end
    if existM1
        ph = iterapp('mldivide',m1fun,m1type,m1fcnstr,p,varargin{:});
        if ~all(isfinite(ph))
            flag = 2;
            resvec = resvec(1:2*ii-1);
            break
        end
    else
        ph = p;
    end
    if existM2
        ph = iterapp('mldivide',m2fun,m2type,m2fcnstr,ph,varargin{:});
        if ~all(isfinite(ph))
            flag = 2;
            resvec = resvec(1:2*ii-1);
            break
        end
    end
    v = iterapp('mtimes',afun,atype,afcnstr,ph,varargin{:});
    rtv = zeros(size(rho)); 
    for i=1:ns
    rtv(i) = rt(:,i)' * v(:,i);
    end
    if (norm(rtv) == 0) | isinf(rtv)
        flag = 4;
        resvec = resvec(1:2*ii-1);
        break
    end
    
    alpha = rho./rtv;
    if isinf(alpha)
        flag = 4;
        resvec = resvec(1:2*ii-1);
        break
    end
    
    if abs(alpha)*norm(ph) < eps*norm(x)
        stag = stag + 1;
    else
        stag = 0;
    end
    
    xhalf = x + ph * diag(alpha);        % form the "half" iterate
    s = r - v * diag(alpha);             % residual associated with xhalf
    normr = norm(s);
    normr_act = normr;
    resvec(2*ii) = normr;
    
    % check for convergence
    if (normr <= tolb || stag >= maxstagsteps || moresteps)
        s = b - iterapp('mtimes',afun,atype,afcnstr,xhalf,varargin{:});
        normr_act = norm(s);
        resvec(2*ii) = normr_act;
        if normr_act <= tolb
            x = xhalf;
            flag = 0;
            iter = ii - 0.5;
            resvec = resvec(1:2*ii);
            break
        else
            if stag >= maxstagsteps && moresteps == 0
                stag = 0;
            end
            moresteps = moresteps + 1;
            if moresteps >= maxmsteps
                if ~warned
                    warning('MATLAB:bicgstab:tooSmallTolerance', ...
                        strcat('Input tol may be smaller than eps*cond(A)',...
                        ' and might not be achieved by BICGSTAB\n',...
                        '         Try to use a bigger tolerance'));
                end
                flag = 3;
                x = xhalf;
                resvec = resvec(1:2*ii);
                break;
            end
        end
    end
    
    if stag >= maxstagsteps
        flag = 3;
        resvec = resvec(1:2*ii);
        break
    end
    
    if normr_act < normrmin        % update minimal norm quantities
        normrmin = normr_act;
        xmin = xhalf;
        imin = ii - 0.5;
    end
    
    if existM1
        sh = iterapp('mldivide',m1fun,m1type,m1fcnstr,s,varargin{:});
        if ~all(isfinite(sh))
            flag = 2;
            resvec = resvec(1:2*ii);
            break
        end
    else
        sh = s;
    end
    if existM2
        sh = iterapp('mldivide',m2fun,m2type,m2fcnstr,sh,varargin{:});
        if ~all(isfinite(sh))
            flag = 2;
            resvec = resvec(1:2*ii);
            break
        end
    end
    t = iterapp('mtimes',afun,atype,afcnstr,sh,varargin{:});
    
    tt = zeros(size(rho)); omega = tt; 
    for j=1:ns
        tt(j) = t(:,j)' * t(:,j);
    
    if (tt(j) == 0) || isinf(tt(j))
        flag = 4;
        resvec = resvec(1:2*ii);
        break
    end
    
    omega(j) = (t(:,j)' * s(:,j)) / tt(j);
    end
    
    if isinf(omega)
        flag = 4;
        resvec = resvec(1:2*ii);
        break
    end
    
    if max(omega)*norm(sh) < eps*norm(xhalf)
        stag = stag + 1;
    else
        stag = 0;
    end
    
    x = xhalf + sh * diag(omega);        % x = (x + alpha * ph) + omega * sh
    r = s - t * diag(omega);
    normr = norm(r);
    normr_act = normr;
    resvec(2*ii+1) = normr;
    
    % check for convergence        
    if (normr <= tolb || stag >= maxstagsteps || moresteps)
        r = b - iterapp('mtimes',afun,atype,afcnstr,x,varargin{:});
        normr_act = norm(r);
        resvec(2*ii+1) = normr_act;
        if normr_act <= tolb
            flag = 0;
            iter = ii;
            resvec = resvec(1:2*ii+1);
            break
        else
            if stag >= maxstagsteps && moresteps == 0
                stag = 0;
            end
            moresteps = moresteps + 1;
            if moresteps >= maxmsteps
                if ~warned
                    warning('MATLAB:bicgstab:tooSmallTolerance', ...
                        'Input tol may be smaller than eps*cond(A) and might not be achieved by BICGSTAB\n         Try to use a bigger tolerance');
                end
                flag = 3;
                resvec = resvec(1:2*ii+1);
                break;
            end
        end        
    end
    
    if normr_act < normrmin        % update minimal norm quantities
        normrmin = normr_act;
        xmin = x;
        imin = ii;
    end
    
    if stag >= maxstagsteps
        flag = 3;
        resvec = resvec(1:2*ii+1);
        break
    end
    
end                                % for ii = 1 : maxit

% returned solution is first with minimal residual
if flag == 0
    relres = normr_act / n2b;
else
    r = b - iterapp('mtimes',afun,atype,afcnstr,xmin,varargin{:});
    if norm(r) <= normr_act
        x = xmin;
        iter = imin;
        relres = norm(r) / n2b;
    else
        iter = ii;
        relres = normr_act / n2b;
    end
end

% only display a message if the output flag is not used
if nargout < 2
    itermsg('bicgstab',tol,maxit,ii,flag,iter,relres);
end