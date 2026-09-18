function KHSS = HSS1D_fsym_to_nsym(KHSS)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        KHSS = HSS1D_fsym_to_nsym(KHSS)
%
%    DESCRIPTION:
%    This function augments the full-symmetric HSS1D input by matching incoming and outgoing
%    skeletons, hence allowing nsym treatment.
%
%    INPUT:
%        KHSS is an HSS structure of a full-symmetric matrix A.
%
%    OUTPUT:
%        KHSS is a "non-symmetric" HSS structure of A.
%


% Match incoming and outgoing skeletons
KHSS(22,:) = KHSS(20,:); % I^{sk} 
KHSS(32,:) = KHSS(30,:); % T
KHSS(33,:) = KHSS(31,:); % J
