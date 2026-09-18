function [J,dr_weights,iistencil] = DR_stencil(h,rk,n,order)
%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%


if (order ~= 0)
  %%% Construct the vector of Duan-Rokhlin weights.
  DVEC = LOCAL_compute_duan_rokhlin_weights(h,rk);
  size(DVEC)
  size([1;1/(h*h);1/(h^4);1/(h^4);1/(h^6);1/(h^6)])
  DVEC_scaled = DVEC.*[1;1/(h*h);1/(h^4);1/(h^4);1/(h^6);1/(h^6)];
end

if (order == 0)
  iistencil  = [0;0;1];
  J          = 0;
  dr_weights = 0;
elseif (order == 4)
  iistencil  = [0;0;1];
  J          = 0;
  dr_weights = DVEC(1);
elseif (order == 6)
   iistencil = [0,  0,  1,  0, -1;...
                0, -1,  0,  1,  0;...
                1,  2,  2,  2,  2];
   J         = iistencil(2,:) + n*iistencil(1,:);
  DR6        = [DVEC_scaled(1) - 2*DVEC_scaled(2); (1/2)*DVEC_scaled(2)];
  dr_weights = DR6(iistencil(3,:)).';
elseif (order == 8)
   iistencil = [0,  0,  1,  0, -1,  0,  2,  0, -2,  1,  1, -1, -1;...
                0, -1,  0,  1,  0, -2,  0,  2,  0  -1,  1,  1, -1;...
                1,  2,  2,  2,  2,  3,  3,  3,  3,  4,  4,  4,  4];
   J         = iistencil(2,:) + n*iistencil(1,:);
   DR8       = [1,  -5/2,  1/2,    1;...
               0,   2/3, -1/6, -1/2;...
               0, -1/24, 1/24,    0;...
               0,     0,    0,  1/4]*DVEC_scaled(1:4);
  dr_weights = DR8(iistencil(3,:)).';
elseif (order == 10)
  iistencil = [0,  0,  1,  0, -1,  0,  2,  0, -2,  0,  3,  0, -3,  1,  1, -1, -1,  1,  2,  2,  1, -1, -2, -2, -1;...
               0, -1,  0,  1,  0, -2,  0,  2,  0  -3,  0,  3,  0, -1,  1,  1, -1, -2, -1,  1,  2,  2,  1, -1, -2;...
               1,  2,  2,  2,  2,  3,  3,  3,  3,  4,  4,  4,  4,  5,  5,  5,  5,  6,  6,  6,  6,  6,  6,  6,  6];
  J         = iistencil(2,:) + n*iistencil(1,:);
  DR10 = [1, -49/18,    7/9,    3/2,  -1/18,  -1/2;...
          0,    3/4, -13/48, -19/24,   1/48,  7/24;...
          0,  -3/40,   1/12,   1/24, -1/120, -1/24;...
          0,  1/180, -1/144,      0,  1/720,     0;...
          0,      0,      0,   5/12,      0,  -1/6;...
          0,      0,      0,  -1/48,      0,  1/48]*DVEC_scaled;
  dr_weights = DR10(iistencil(3,:)).';
else
  fprintf(1,'ERROR: This option for "order" is not implemented.\n')
  %keyboard
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DVEC = LOCAL_compute_duan_rokhlin_weights(h,rk)

fid = fopen('indata.txt','w');
fprintf(fid,'h       = %23.15e\n',h);
fprintf(fid,'real_rk = %23.15e\n',real(rk));
fprintf(fid,'imag_rk = %23.15e\n',imag(rk));

fclose(fid);

wpath = strsplit(pwd,'HSSC_DirectSolver');                      % (1) DEFAULT PATH/EXEC
wpath = strcat(char(wpath(1)),'HSSC_DirectSolver/Common/');
if strcmp(computer,'MACI64')
    wpath = strcat(wpath,'FORTRAN_get_DR_weights_MAC');  
elseif strcmp(computer,'PCWIN') || strcmp(computer,'PCWIN64')
    wpath = strcat(wpath,'FORTRAN_get_DR_weights_WINDOWS.exe');
elseif strcmp(computer,'GLNXA64')
    wpath = strcat(wpath,'FORTRAN_get_DR_weights_LINUX');
end
% wpath = 'path/execuatable'                                    % (2) SELF-COMPILED EXE
system(wpath);                                                  % (3) RUN EXECUTABLE

fid = fopen('outdata.txt');

mydata = textscan(fid,'%s %s %d',1);
ier = mydata{3};

mydata = textscan(fid,'%s %s %f',1);
h = mydata{3};

mydata = textscan(fid,'%s %s %f %f',1);
real_rk = mydata{3};
imag_rk = mydata{4};

if (~(ier == 0) | (abs(rk - real_rk - 1i*imag_rk) > 1e-8))
  fprintf(1,['ERROR: the external function "compute_duan_rokhlin_weights" ',...
             'did not complete correctly.\n'])
  %keyboard
end

mydata = textscan(fid,'%s %s %f %f',6);
DVEC = mydata{3} + 1i*mydata{4};

fclose(fid);
delete('indata.txt','outdata.txt')    

return