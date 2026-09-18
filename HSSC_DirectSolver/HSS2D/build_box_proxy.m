function [X_proxy,X_extra] = build_box_proxy(cB,colrad,rowrad,layers,h)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%     FUNCTION CALL:
%         [X_proxy,X_extra] = build_box_proxy(cB,colrad,rowrad,layers,h)
%
%     DESCRIPTION:
%         This function constructs the proxy of exterior points surrounding a box. By placing
%         point charges (discrete analog of an equivalent density) around the boundary, we can
%         represent the interaction between the target points (inside the box) and the source
%         points (outside the box). This allows for the accelerated compression of off-diagonal
%         block rows and columns of the HSS matrix.
%
%         The width of the layer of proxy points depends on the accuracy requested. We found for
%         the Laplace kernel, a layer of width 1 and width 2 leads to relative accuracy of about
%         10^-5 and 10^-10, respectively. For the Helmholtz kernel similar accuracy was
%         observed, however, thicker skeleton layers are recommended to avoid problems
%         associated with resonances.
%
%
%                     * * * * * * * * * *         * proxy points (eg. 2 layers)
%                     * * * * * * * * * *
%                     * *             * *
%                     * *             * *
%                     * *             * *
%                     * *             * *
%                     * *             * *
%                     * * * * * * * * * *
%                     * * * * * * * * * *
%
%     INPUT:
%         cB      <1x2 float>     Center of the box.
%         colrad  <float>         Half the width of the box.
%         rowrad  <float>         Half the height of the box.
%         layers  <int>           Number of layers around the box.
%         h       <float>         Separation between proxy points (both dimensions).
%
%
%     OUTPUT:
%         X_proxy <mx2 float>         Points in the proxy.
%         X_extra <(m-k)x2 float>     Corner points are removed to make the proxy the same
%                                         size as the skeleton set (see Remark 3.4 of the paper).
%




width = (layers-1/2)*h; 
[Xp1,Yp1] = meshgrid(cB(1) - colrad - width:h:cB(1) - colrad - h/2 , cB(2) - rowrad - width:h:cB(2) + rowrad + width); 
[Xp2,Yp2] = meshgrid(cB(1) - colrad + h/2:h:cB(1) + colrad - h/2 , cB(2) + rowrad + width:-h:cB(2) + rowrad + h/2); 
[Xp3,Yp3] = meshgrid(cB(1) + colrad + h/2:h:cB(1) + colrad + width , cB(2) + rowrad + width:-h: cB(2) - rowrad - width);
[Xp4,Yp4] = meshgrid(cB(1) + colrad - h/2:-h:cB(1) - colrad + h/2 , cB(2) - rowrad - h/2:-h:cB(2) - rowrad - width); 
Xp1 = Xp1'; Yp1 = Yp1'; 
Xp3 = Xp3'; Yp3 = Yp3'; 

X1 = [Xp1(:) Yp1(:)]; X2 = [Xp2(:) Yp2(:)]; X3 = [Xp3(:) Yp3(:)]; X4 = [Xp4(:) Yp4(:)];
clear Xp1 Yp1 Xp2 Yp2 Xp3 Yp3 Xp4 Yp4; 
X_proxy = [X1 ; X2 ; X3 ; X4];

M = size(X_proxy,1); 

% Option to remove points from proxy to make it of the same size as the
% skeleton set. 

%if M>250
    % Points to remove per corner
    r1 = layers^2-0.5; 
    r2 = layers^2+0.5;  
    
    %c = (M-m)/4; 
    %r1 = (c-1)/2; 
    %r2 = (c+1)/2;
    
    % Indices to be removed
    I1 = (abs(X1(:,1)-cB(1))>(colrad + width - h/2) & (X1(:,2)-cB(2)>(rowrad + width-r1*h) | X1(:,2)-cB(2)<-(rowrad + width-r2*h))) | abs(X1(:,2)-cB(2))>(rowrad + width - h/2);
    I2 = (abs(X2(:,2)-cB(2))>(rowrad + width - h/2) & (X2(:,1)-cB(1)>(colrad + width-r1*h) | X2(:,1)-cB(1)<-(colrad + width-r2*h)));
    I3 = (abs(X3(:,1)-cB(1))>(colrad + width - h/2) & (X3(:,2)-cB(2)>(rowrad + width-r2*h) | X3(:,2)-cB(2)<-(rowrad + width-r1*h))) | abs(X3(:,2)-cB(2))>(rowrad + width - h/2);
    I4 = (abs(X4(:,2)-cB(2))>(rowrad + width - h/2) & (X4(:,1)-cB(1)>(colrad + width-r2*h) | X4(:,1)-cB(1)<-(colrad + width-r1*h)));
    
    X_proxy = [X1(~I1,:) ; X2(~I2,:) ; X3(~I3,:) ; X4(~I4,:)]; 
    X_extra = [X1(I1,:)  ; X2(I2,:)  ; X3(I3,:)  ; X4(I4,:)]; 
%else
%    X_extra = zeros(size(X_proxy,1),0);
%end


