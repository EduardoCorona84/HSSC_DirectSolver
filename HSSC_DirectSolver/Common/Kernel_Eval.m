function A = Kernel_Eval(X1,X2,params)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%     FUNCTION CALL:
%         A = Kernel_Eval(X1,X2,params)
%
%     DESCRIPTION:
%         This function takes 2D point arrays X1 and X2, and returns the matrix
%         A = K[X1,X2]. The only part of params used is params.flagpot, which is a name for the
%         kernel function being to be used. A few default kernels are available but you can also
%         input your own kernel function (see RUNNING THE CODE section).
%
%     INPUT:
%         X1      <mx2 float>     Source points
%         X2      <mx2 float>     Target points (size need not be the same as X1)
%
%         params.flagpot -    'SL_H_2D' - 2D Helmholtz
%                             'SL_H_3D' - 3D Helmholtz
%                             'SL_L_2D' - 2D Laplace
%                             'SL_L_3D' - 3D Laplace
%                             'SL_Y_2D' - 2D Yukawa
%                             'SL_Y_3D' - 3D Yukawa
%

h = params.h; 

% X and Y grid
[Y_g1,  X_g1  ] = meshgrid(X2(:,1), X1(:,1));
[Y_g2,  X_g2  ] = meshgrid(X2(:,2), X1(:,2));

% Determine functions b(x) and c(y)
[b,c] = LOCAL_get_bc(X_g1,X_g2,Y_g1,Y_g2,params,4);

% 2D (plane or curves on plane)
 if params.dim == 2
     % den = ||X-Y||^2
     den = (X_g1 - Y_g1).^2 + (X_g2 - Y_g2).^2;
     
     % 2D Helmholtz
    if strcmp(params.flag_pot,'SL_H_2D')
        kh = params.kh; 
        C = besselh(0,kh);
        if params.order == 4
            w = 1 + kh^2*(0.25*1i)*params.dr_weights; 
            A = w*(den==0) + (h^2)*kh^2*(0.25*1i)*b.*(besselh(0,kh*sqrt(den + (den==0))) - C*(den==0)).*c;
        else
            % Have to add q-diagonal correction here. 
            A = (den==0) + (h^2)*kh^2*(0.25*1i)*b.*(besselh(0,kh*sqrt(den + (den==0))) - C*(den==0)).*c;
        end
    % 3D Helmholtz    
    elseif strcmp(params.flag_pot,'SL_H_3D')
        kh = params.kh; 
        C = exp(1i*kh);
        A = (den==0) + (h^2)*kh*(exp(1i*kh*sqrt(den + (den==0)))./sqrt(den + (den==0)) - C*(den==0)); 
    % 2D Laplace    
    elseif strcmp(params.flag_pot,'SL_L_2D')
        A = (den==0) + (h^2)*0.5*b.*log((den + (den==0))).*c;
    % 3D Laplace    
    elseif strcmp(params.flag_pot,'SL_L_3D')
        A = (den==0) + (h^2)*b.*(1./sqrt(den + (den==0)) - (den==0)).*c; 
    % 2D Yukawa    
    elseif strcmp(params.flag_pot,'SL_Y_2D')
        C = besselk(0,params.kh);  
        A = (den==0) + (h^2)*(0.5/pi)*b.*(besselk(0,params.kh*(den + (den==0))) - C*(den==0)).*c; 
    %3D Yukawa    
    elseif strcmp(params.flag_pot,'SL_Y_3D')
        A = (den==0) + (h^2)*b.*(exp(-params.kh*den)./sqrt(den + (den==0)) - (den==0)).*c;  
    elseif strcmp(params.flag_pot,'Gauss_2D')
        s = 4*h^2; 
        A = (1/(4*pi*h^2))*exp(-(1/s)*den); 
    elseif strcmp(params.flag_pot,'Kernel_Name')
            % define distance function den
            % define kernel A = (h^2)/2 b(X) K(|X-Y|) C(Y)
        % den = (X_g1-Y_g1).^2 + (X_g2-Y_g2).^2;          
        % A = (h^2)/2 * b(X_g1,X_g2) K(den) c(Y_g1,Y_g2);
    end
 % 3D (surfaces and volume)   
 else
     % Z grid
     Z1 = X1(:,3); Z2 = X2(:,3); 
     [Y_g3,X_g3] = meshgrid(Z2,Z1);
     % den = ||X-Y||^2
     den = (X_g1 - Y_g1).^2 + (X_g2 - Y_g2).^2 + (X_g3 - Y_g3).^2;
     
     if strcmp(params.flag_pot,'SL_L_3D')  
        A = (den==0) + (h^2)*(1./sqrt(den + (den==0)) - (den==0)); 
     elseif strcmp(params.flag_pot,'Gauss_3D')
        s = 4*h^2; 
        A = (1/(4*pi*h^2))*exp(-(1/s)*den);
     end
 end
end

 % b function for NTI Laplace example
function b = b_func(xx1,xx2,x0,y0,alpha)

if alpha==0
    b = ones(size(xx1)); 
else
    dd1 = (xx1 - x0);
    dd2 = (xx2 - y0);
    b   = 1 + alpha*exp( - dd1.*dd1 - dd2.*dd2);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [b,c] = LOCAL_get_bc(X_g1,X_g2,Y_g1,Y_g2,params,opts)
   
h = params.h; 

if params.sym == 0 && params.kh>0
   % Option for non-symmetric Lippmann-Schwinger 
   if params.proxy <= 0
        b = Bump_function(X_g1,X_g2,opts); 
   else
        b = ones(size(X_g1)); 
   end
      
   if params.proxy >= 0 
        c = ones(size(Y_g1)); 
        %c = Bump_function(Y_g1,Y_g2,opts);
   else
        c = ones(size(Y_g1)); 
   end
elseif params.transinv == 0 
   % Option for symmetric Laplace / Lippmann-Schwinger 
   if params.proxy <= 0
        if params.kh>0
             b = sqrt(Bump_function(X_g1,X_g2,opts,h)); 
        else
             b = b_func(X_g1,X_g2,0.3,0.6,0.25); 
        end
   else
        b = ones(size(X_g1)); 
   end
       
   if params.proxy >=0
        if params.kh>0
            c = sqrt(Bump_function(Y_g1,Y_g2,opts,h));
        else
            if params.sym
                c = b_func(Y_g1,Y_g2,0.3,0.6,0.25);
            else
                c = b_func(Y_g1,Y_g2,0.2,-0.6,0.25);
            end
        end
   else
        c = ones(size(Y_g1)); 
   end  
else
   % Else, b and c are ones
   b = ones(size(X_g1)); 
   c = ones(size(Y_g1)); 
end
end