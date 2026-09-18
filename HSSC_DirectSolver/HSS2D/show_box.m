function show_box(cB,lev,color)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%     FUNCTION CALL:
%         show_box(cB,lev,color)
%
%     DESCRIPTION:
%         This function allows us to plot a box with center cB at levels lev = [lev_x lev_y]
%         (where level zero corresponds to radius 1).
%




linx = [cB(1)-2^(-lev(1)) , cB(1)+2^(-lev(1))];
liny = [cB(2)-2^(-lev(2)) , cB(2)+2^(-lev(2))];

hold on
plot((cB(1)-2^(-lev(1)))*ones(1,2),liny,color,'LineWidth',4)
plot((cB(1)+2^(-lev(1)))*ones(1,2),liny,color,'LineWidth',4)
plot(linx,(cB(2)-2^(-lev(2)))*ones(1,2),color,'LineWidth',4)
plot(linx,(cB(2)+2^(-lev(2)))*ones(1,2),color,'LineWidth',4)


