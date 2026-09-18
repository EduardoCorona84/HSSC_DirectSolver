function KNEW = HSS1D_nsym_to_fsym(KHSS,acc)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        KNEW = HSS1D_nsym_to_fsym(KHSS,acc)
%
%    DESCRIPTION:
%        This function converts the non-symmetric HSS matrix KHSS into a full-symmetric HSS matrix
%        KNEW. This process corresponds to computing (A+A')/2 with HSS1D sum.
%
%    INPUT:
%        KHSS is the {50 x nboxes} cell array that encodes the HSS matrix to be converted.
%        acc (double) is the desired accuracy.
%
%    OUTPUT:
%        KNEW is the {50 x nboxes} cell array of the converted HSS matrix.
%


global BOX 

nboxes = size(KHSS,2);
KNEW = KHSS; 

% Compress all KHSS, going from smaller to larger.
for ibox = nboxes:(-1):1
    
    if ibox>1
    k0 = KHSS{BOX.KSKEL,ibox};  
    R = [eye(k0)  KHSS{BOX.T_UP,ibox}]; R(:,KHSS{BOX.J_UP,ibox}) = R;
    Lt = [eye(k0)  KHSS{BOX.T_DN,ibox}]; Lt(:,KHSS{BOX.J_DN,ibox}) = Lt; 
    
  if ( (KHSS{BOX.C1,ibox}<=0) & (KHSS{BOX.C2,ibox}<=0) ) % ibox has no sons.
    ind_skel_out = KHSS{BOX.END1,ibox} - 1 + (1:KHSS{BOX.END2,ibox});
  elseif ( (KHSS{BOX.C1,ibox} > 0) & (KHSS{BOX.C2,ibox}<= 0) ) % ibox has a left son
    ison1        = KHSS{BOX.C1,ibox}; k1 = KHSS{BOX.KSKEL,ison1};
    ind_skel_out = KNEW{BOX.I_SKUP,ison1};
    R = R*KHSS{BOX.RHO,ison1}(1:k1,:);
    Lt = Lt*KHSS{BOX.RHO,ison1}(k1+1:end,:);
  elseif ( (KHSS{BOX.C1,ibox}<=0) & (KHSS{BOX.C2,ibox}>0) ) % ibox has a right son
    ison2        = KHSS{BOX.C2,ibox}; k2 = KHSS{BOX.KSKEL,ison2}; 
    ind_skel_out = KNEW{BOX.I_SKUP,ison2};
    R = R*KHSS{BOX.RHO,ison2}(1:k2,:);
    Lt = Lt*KHSS{BOX.RHO,ison2}(k2+1:end,:);
  else % ibox has two sons
    ison1        = KHSS{BOX.C1,ibox}; k1 = KHSS{BOX.KSKEL,ison1}; 
    ison2        = KHSS{BOX.C2,ibox}; k2 = KHSS{BOX.KSKEL,ison2}; 
    ind_skel_out = [KNEW{BOX.I_SKUP,ison1},KNEW{BOX.I_SKUP,ison2}];
    
    R = R*[KHSS{BOX.RHO,ison1}(1:k1,:) zeros(k1,KNEW{BOX.KSKEL,ison2}) ; zeros(k2,KNEW{BOX.KSKEL,ison1}) KHSS{BOX.RHO,ison2}(1:k2,:)];
    Lt = Lt*[KHSS{BOX.RHO,ison1}(k1+1:end,:) zeros(k1,KNEW{BOX.KSKEL,ison2}) ; zeros(k2,KNEW{BOX.KSKEL,ison1}) KHSS{BOX.RHO,ison2}(k2+1:end,:)];
  end
  
  RC = [R ; Lt];
  
  % Compute the skeletons.
  [T,J] = ID(RC,acc);
  k = size(T,1);
 
  % Record the outgoing skeletons:
  KNEW{BOX.I_SKUP,ibox} = ind_skel_out(J(1:k));
  KNEW{BOX.T_UP,ibox} = T;
  KNEW{BOX.J_UP,ibox} = J;
  % Record the incoming skeletons:
  KNEW{BOX.I_SKDN,ibox} = [];
  KNEW{BOX.T_DN,ibox} = [];
  KNEW{BOX.J_DN,ibox} = [];
  % Record nskel and kskel
  KNEW{BOX.NSKEL,ibox} = length(ind_skel_out);
  KNEW{BOX.KSKEL,ibox} = k;
  KHSS{BOX.RHO,ibox} = RC(:,J(1:k));
    end
  
  if ( (KHSS{BOX.C1,ibox}>0) & (KHSS{BOX.C2,ibox}>0) )
      ison1        = KHSS{BOX.C1,ibox}; k1 = KHSS{BOX.KSKEL,ison1}; 
      ison2        = KHSS{BOX.C2,ibox}; k2 = KHSS{BOX.KSKEL,ison2}; 
      KNEW{BOX.M_SIB,ison1} = (KHSS{BOX.RHO,ison1}(k1+1:end,:).'*KHSS{BOX.M_SIB,ison1}*KHSS{BOX.RHO,ison2}(1:k2,:) +...
                            KHSS{BOX.RHO,ison1}(1:k1,:).'*KHSS{BOX.M_SIB,ison2}.'*KHSS{BOX.RHO,ison2}(k2+1:end,:))/2;
      KNEW{BOX.M_SIB,ison2} = KNEW{BOX.M_SIB,ison1}.';
  elseif ( (KNEW{BOX.C1,ibox}<=0) & (KNEW{BOX.C2,ibox}<=0) )
      KNEW{BOX.M_SELF,ibox} = (KHSS{BOX.M_SELF,ibox} + KHSS{BOX.M_SELF,ibox}.')/2;
  end
  
end