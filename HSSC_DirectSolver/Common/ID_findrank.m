function [krank,RA] = ID_findrank(m,n,matvecA,dist,k0,eps)
%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%   DESCRIPTION:
%      Estimates the numerical rank krank of a matrix a to precision
%      eps, where the routine matvect applies the transpose of a
%      to an arbitrary vector. This routine applies the transpose of a
%      to krank random vectors, and returns the resulting vectors
%      as the columns of ra.  
% 
%   INPUT:
%         eps -- precision defining the numerical rank
%         m -- first dimension of a
%         n -- second dimension of a
%         matvectA -- routine which applies the transpose
%                    of the matrix whose rank is to be estimated
%                    to an arbitrary vector; this routine must have
%                    a calling sequence of the form
% 
%                    matvect(m,x,n,y,p1,p2,p3,p4),
% 
%                    where m is the length of x,
%                    x is the vector to which the transpose
%                    of the matrix is to be applied,
%                    n is the length of y,
%                    y is the product of the transposed matrix and x,
%                    and p1, p2, p3, and p4 are user-specified parameters
%         p1 -- parameter to be passed to routine matvect
%         p2 -- parameter to be passed to routine matvect
%         p3 -- parameter to be passed to routine matvect
%         p4 -- parameter to be passed to routine matvect
% 
%     OUTPUT:
%         krank -- estimate of the numerical rank of a
%         ra -- product of the transpose of a and a matrix whose entries
%               are pseudorandom realizations of i.i.d. random numbers,
%               uniformly distributed on [0,1];
%               ra must be at least 2*n*krank real*8 elements long
%


krank = 0; 
residual = 1; 
scal = zeros(n,1); 
num = k0; 

while (residual>eps && krank<min(m,n))
    if strcmp(dist,'unif')
        x = rand(m,num);
    elseif strcmp(dist,'norm')
        x = randn(m,num);
    else
        x = ones(m,num); 
    end
   
   y(:,krank+(1:num)) = matvecA(x); 
   RA(:,krank+(1:num)) = y(:,krank+(1:num)); 
   
   i=1;
   while i<=num
   if krank>0
       
       ifrescal = 0; 
       for k=0:krank-1
          [y(k+1:n,krank+1),scal(k+1)] = ID_houseapp(n-k,vn(1:n-k,k+1),y(k+1:n,krank+1),ifrescal,scal(k+1)); 
       end
       
   end
   
       [residual,vn(1:n-krank,krank+1),scal(krank+1)] = ID_house(n-krank,y(krank+1:n,krank+1));
       residual = abs(residual); 
       
       krank = krank + 1;
       if ~(residual>eps && krank<min(m,n))
       i = num+1; 
       else
           i=i+1;
       end   
   end
   num = 10;
end

krank = krank-1; 