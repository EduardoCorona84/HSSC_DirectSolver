function [FF,xxsupp] = Bump_function(XX1,XX2,flag_option,h)
%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%


if (nargin == 2)
  flag_option = 4;
end

if nargin == 3
   h = 0; 
end

xx1 = XX1(:);
xx2 = XX2(:);

if (flag_option == 1)
  %%% Tensor product of exp(-1/(1-x^2)) bumps.

  ii = find( (abs(xx1)<1) & (abs(xx2)<1));
  ff = zeros(size(xx1));
  ff(ii) = exp(-1./(1 - xx1(ii).*xx1(ii))).*...
           exp(-1./(1 - xx2(ii).*xx2(ii)));
  FF = reshape(ff,size(XX1));
  xxsupp = [-1,1,1,-1,-1;-1,-1,1,1,-1];
     
elseif (flag_option == 2)
  %%% Rotated bump exp(-1/(1 - r^2))
    
  A = [1,0.4;...
       0,1];
  yy1 = A(1,1)*xx1 + A(1,2)*xx2;
  yy2 = A(2,1)*xx1 + A(2,2)*xx2;
  rrsq = yy1.*yy1 + yy2.*yy2;
  ii  = find( rrsq < 1 );
  ff  = zeros(size(xx1));
  ff(ii) = exp(-1./(1 - rrsq(ii)));
  FF = reshape(ff,size(XX1));
  tt = linspace(0,2*pi,1000);
  xxsupp = inv(A)*[cos(tt);sin(tt)];
  
elseif (flag_option == 3)
  %%% Gaussian
  
  rrsq = xx1.*xx1 + xx2.*xx2;
  ff   = exp(-15*rrsq);
  FF   = reshape(ff,size(XX1));
  xxsupp = NaN*ones(2,5);
elseif (flag_option == 4)
    %%% 1-0 exponential decay w/ tanh
    c = 4*h; %1/14; %4*params.h;     
    delta = 1/10; 
    
    D = 1/delta; 
    ff1 = 0.5*tanh(D*(-abs(xx1)+1-c))+0.5; 
    ff2 = 0.5*tanh(D*(-abs(xx2)+1-c))+0.5; 
    ff = ff1.*ff2; 
    FF = reshape(ff,size(XX1)); 
    xxsupp = NaN*ones(2,5); 
  
end