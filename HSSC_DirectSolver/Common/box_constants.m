%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%	DESCRIPTION:
%		This file defines global constants used for HSS structure indices. It specifies how an 
%		HSS1D matrix is encoded in a cell-array.
%	
%		BOX is a global variable throughout the workspace that defines the indices of an HSS1D 
%		structure consistent with the notation of the paper.
%
%		ibox is the index of a box in the HSS tree.
%		HSS is a cell array {50,nboxes} that encodes the HSS1D form of a matrix as follows.
%			HSS{BOX.LEVEL,ibox} - level
%			HSS{BOX.PARENT,ibox} - parent
%			HSS{BOX.C1:C2,ibox} - children c1 and c2
%			HSS{BOX.END1:END2,ibox} - index endpoints
%			HSS{BOX.NSKEL,ibox} - nskel 
%			HSS{BOX.KSKEL,ibox} - kskel
%			HSS{BOX.I_SKUP,ibox} - I_skup outgoing index with skeleton points
%			HSS{BOX.I_SKDN,ibox} - I_skdn incoming index with skeleton points 
%			HSS{BOX.T_UP,ibox} - T_up outgoing interpolation matrix
%			HSS{BOX.J_UP,ibox} - J_up outgoing skeleton permutation vector
%			HSS{BOX.T_DN,ibox} - T_dn incoming interpolation matrix
%			HSS{BOX.J_DN,ibox} - J_dn incoming skeleton permutation vector
%			HSS{BOX.M_SELF,ibox} - Self-Interaction matrices (leaves only) 
%			HSS{BOX.M_SIB,ibox} - Sibling Interaction matrices (except top)
%

% CONSTANTS FOR CELL ARRAY INDICES
global BOX 
BOX.DATA = 1:7;
BOX.LEVEL = 2;
BOX.PARENT = 3;
BOX.C1 = 4;
BOX.C2 = 5;
BOX.END1 = 6;
BOX.END2 = 7;
BOX.NSKEL = 8;
BOX.KSKEL = 9;
BOX.I_SKUP = 20;
BOX.I_SKDN = 22;
BOX.T_UP = 30;
BOX.J_UP = 31;
BOX.T_DN = 32;
BOX.J_DN = 33;
BOX.RHO = 34;
BOX.LAMBDA = 35;
BOX.M_SELF = 40;
BOX.M_SIB = 46;
