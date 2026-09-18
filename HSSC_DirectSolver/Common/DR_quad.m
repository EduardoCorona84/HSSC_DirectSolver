function [xx,ww] = DR_quad(a,n,xxt,rk,order)
%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%


h = 2*a/(n-1);

%%% Create the actual grid.
[XX2,XX1] = meshgrid(linspace(-a,a,n));
xx = [reshape(XX1,1,numel(XX1));...
      reshape(XX2,1,numel(XX2))];
clear XX1 XX2
  
%%% Verify that xxt is an interior grid point.
[mindist,i] = min(abs(xx(1,:) - xxt(1)) + abs(xx(2,:) - xxt(2)));
if ~(numel(i) == 1)
  fprintf(1,'ERROR: Did not find precisely one min point.\n')
  %keyboard
end
if (mindist > 1e-13)
  fprintf(1,'ERROR: It seems xxt is not a grid point.\n')
  %keyboard
end
if (max(abs(xxt)) > (a - 2.5*h))
  fprintf(1,'ERROR: It seems xxt is too close to the boundary.\n')
  %keyboard
end
    
%%% Set up the vector of "regular weights".
g       = [0.5*h,h*ones(1,n-2),0.5*h];
GG      = g' * g;
gg      = reshape(GG,1,numel(GG));
rrsq    = (xx(1,:)-xxt(1)).^2 + (xx(2,:) - xxt(2)).^2;
rrsq(i) = 1;
ww      = gg.*besselh(0,rk*sqrt(rrsq));

%%% Construct the vector of Duan-Rokhlin weights.
DVEC = LOCAL_compute_duan_rokhlin_weights(h,rk);
DVEC_scaled = DVEC.*[1;1/(h*h);1/(h^4);1/(h^4);1/(h^6);1/(h^6)];

if (order == 0)
  ww(i) = 0;
  return    
elseif (order == 4)
  ww(i) = DVEC(1);
  return
elseif (order == 6)
  I1     = i + [-n,n,-1,1];
  DR6    = [DVEC_scaled(1) - 2*DVEC_scaled(2); (1/2)*DVEC_scaled(2)];
  ww(i)  = DR6(1);
  ww(I1) = ww(I1) + DR6(2)*ones(size(ww(I1)));
  return
elseif (order == 8)
  I1     = i + [-n,n,-1,1];
  I2     = i + [-2*n,2*n,-2,2];
  I3     = i + [-n-1,-n+1,n-1,n+1];
  DR8    = [1,  -5/2,  1/2,    1;...
            0,   2/3, -1/6, -1/2;...
            0, -1/24, 1/24,    0;...
            0,     0,    0,  1/4]*DVEC_scaled(1:4);
  ww(i)  = DR8(1);
  ww(I1) = ww(I1) + DR8(2)*ones(size(ww(I1)));
  ww(I2) = ww(I2) + DR8(3)*ones(size(ww(I2)));
  ww(I3) = ww(I3) + DR8(4)*ones(size(ww(I3)));
  return
elseif (order == 10)
  I1   = i + [-n,n,-1,1];
  I2   = i + [-2*n,2*n,-2,2];
  I3   = i + [-3*n,3*n,-3,3];
  I4   = i + [-n-1,-n+1,n-1,n+1];
  I5   = i + [-2*n-1,-2*n+1,-n-2,-n+2,n-2,n+2,2*n-1,2*n+1];
  DR10 = [1, -49/18,    7/9,    3/2,  -1/18,  -1/2;...
          0,    3/4, -13/48, -19/24,   1/48,  7/24;...
          0,  -3/40,   1/12,   1/24, -1/120, -1/24;...
          0,  1/180, -1/144,      0,  1/720,     0;...
          0,      0,      0,   5/12,      0,  -1/6;...
          0,      0,      0,  -1/48,      0,  1/48]*DVEC_scaled;
  ww(i)  = DR10(1);
  ww(I1) = ww(I1) + DR10(2)*ones(size(ww(I1)));
  ww(I2) = ww(I2) + DR10(3)*ones(size(ww(I2)));
  ww(I3) = ww(I3) + DR10(4)*ones(size(ww(I3)));
  ww(I4) = ww(I4) + DR10(5)*ones(size(ww(I4)));
  ww(I5) = ww(I5) + DR10(6)*ones(size(ww(I5)));
  return
else
  fprintf(1,'ERROR: This option for "order" is not implemented.\n')
  %keyboard
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DVEC = LOCAL_compute_duan_rokhlin_weights(h,rk)

fid = fopen('indata.txt','w');
fprintf(fid,'h       = %23.15e\n',h);
fprintf(fid,'real_rk = %23.15e\n',real(rk));
fprintf(fid,'imag_rk = %23.15e\n',imag(rk));

fclose(fid);

wpath = strsplit(pwd,'HSSC_DirectSolver');                  
wpath = strcat(char(wpath(1)),'HSSC_DirectSolver/Common/');     % (2) PATH OF COMMON FOLDER
if strcmp(computer,'MACI64')                                    % (3) ADD SUFFIX
    wpath = strcat(wpath,'FORTRAN_get_DR_weights_MAC');  
elseif strcmp(computer,'PCWIN') || strcmp(computer,'PCWIN64')
    wpath = strcat(wpath,'FORTRAN_get_DR_weights_WINDOWS.exe');
elseif strcmp(computer,'GLNXA64')
    wpath = strcat(wpath,'FORTRAN_get_DR_weights_LINUX');
end
% wpath = strcat(wpath,'<NAME OF EXECUTABLE>');                 % (4) SELF-COMPILED EXE
system(wpath);                                                  % (5) RUN EXECUTABLE

fid = fopen('outdata.txt');

mydata = textscan(fid,'%s %s %d',1);
ier = mydata{3};

mydata = textscan(fid,'%s %s %f',1);
h = mydata{3};

mydata = textscan(fid,'%s %s %f %f',1);
real_rk = mydata{3};
imag_rk = mydata{4};

if (~(ier == 0) || (abs(rk - real_rk - 1i*imag_rk) > 1e-8))
  fprintf(1,['ERROR: the external function "compute_duan_rokhlin_weights" ',...
             'did not complete correctly.\n'])
  %keyboard
end

mydata = textscan(fid,'%s %s %f %f',6);
DVEC = mydata{3} + 1i*mydata{4};

fclose(fid);
delete('indata.txt','outdata.txt')

return