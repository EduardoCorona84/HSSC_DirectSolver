function [krank,Y] = QR_findrank(m,n,matvecA,dist,k0,C,eps)
%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%     estimates the numerical rank krank of a matrix a to precision
%     eps, where the routine matvect applies the transpose of a
%     to an arbitrary vector. This routine applies the transpose of a
%     to krank random vectors, and returns the resulting vectors
%     as the columns of ra.
%     
%     input:
%     eps -- precision defining the numerical rank
%     m -- first dimension of a
%     n -- second dimension of a
%     matvectA -- routine which applies the transpose
%               of the matrix whose rank is to be estimated
%               to an arbitrary vector; this routine must have
%               a calling sequence of the form
%     
%               matvect(m,x,n,y,p1,p2,p3,p4),
%     
%               where m is the length of x,
%               x is the vector to which the transpose
%               of the matrix is to be applied,
%               n is the length of y,
%               y is the product of the transposed matrix and x,
%               and p1, p2, p3, and p4 are user-specified parameters
%     p1 -- parameter to be passed to routine matvect
%     p2 -- parameter to be passed to routine matvect
%     p3 -- parameter to be passed to routine matvect
%     p4 -- parameter to be passed to routine matvect
%     
%     output:
%     krank -- estimate of the numerical rank of a
%     ra -- product of the transpose of a and a matrix whose entries
%          are pseudorandom realizations of i.i.d. random numbers,
%          uniformly distributed on [0,1];
%          ra must be at least 2*n*krank real*8 elements long
%


k0 = min(k0,min(m,n)); 

if strcmp(dist,'unif')
    G = rand(m,k0); 
elseif strcmp(dist,'norm')
    G = randn(m,k0); 
end
Y = matvecA(G); 
[Q,R,~] = qr(Y,0); 

D = abs(diag(R));
if length(D)>0 
residual = D(end);
 else
   residual = 0; 
end
 
krank = sum(D>eps);

while (residual>eps && krank<min(m,n))
    num = min(C,min(m,n)-krank); 
    
    if strcmp(dist,'unif')
        G = rand(m,num);
    elseif strcmp(dist,'norm')
        G = randn(m,num);
    end
   
   Y(:,krank+(1:num)) = matvecA(G); 
   M = Y(:,krank+(1:num)) - Q*(Q'*Y(:,krank+(1:num))); 
   
   [Qc,R,~] = qr(M,0); 
   D = abs(diag(R)); 
   if length(D)>0
   residual = D(end);
   else 
   residual = 0; 
   end
 
   krank = krank + sum(D>eps);
   
   Qc = Qc - Q*(Q'*Qc);
   Q = [Q Qc];
end

%krank = krank-1; 
Y = Y(:,1:krank);
%krank = krank-1; 
