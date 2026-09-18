function uu = HSS1D_apply_nsym(K_HSS,qq,dim)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2008-09 P.G. Martinsson (C) 2011-2013 E. Corona, P.G. Martinsson, D. Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        uu = HSS1D_apply_nsym(K_HSS,qq,dim)
%
%    DESCRIPTION:
%        This function performs a fast matrix-vector multiply of an square, non symmetric HSS1D
%        matrix and a vector : uu=Aqq
%
%    INPUT:
%        A denotes the matrix encoded in HSS1D structure K_HSS.
%        qq is the vector to be left-multiplied by K_HSS.
%        dim (int) is the dimension of the kernal
%
%    OUTPUT:
%        uu = Aqq vector.
%


global BOX 

if nargin == 2
    dim = 1; 
end

nboxes = size(K_HSS,2);
FIELDS = cell(2,nboxes);
Ntot = length(qq)/dim;

% Initialize the fields on the leaves
%  - Construct the outgoing fields.
%  - Set the incoming fields to zero.
for ibox = nboxes:(-1):2
    
  if ( (K_HSS{BOX.C1,ibox}<=0) && (K_HSS{BOX.C2,ibox}<=0) )
    ind            = dim*K_HSS{BOX.END1,ibox} - dim + (1:dim*K_HSS{BOX.END2,ibox});
    FIELDS{1,ibox} = LOCAL_skel_proj(K_HSS{BOX.T_UP,ibox},K_HSS{BOX.J_UP,ibox},qq(ind,:));
    FIELDS{2,ibox} = zeros(length(K_HSS{BOX.I_SKDN,ibox}),size(qq,2));
  end
end

% Construct all outgoing potentials:
for ibox = nboxes:(-1):2
  if ( (K_HSS{BOX.C1,ibox}>0) || (K_HSS{BOX.C2,ibox}>0) ) % ibox has at least one son.
      if (K_HSS{BOX.C1,ibox}>0)
          ison1 = K_HSS{BOX.C1,ibox};
          if (K_HSS{BOX.C2,ibox}>0)
              ison2 = K_HSS{BOX.C2,ibox};
              qloc  = [FIELDS{1,ison1};FIELDS{1,ison2}];
          else
              qloc  = [FIELDS{1,ison1}];
          end
      else
          ison2 = K_HSS{BOX.C2,ibox};
          qloc  = [FIELDS{1,ison2}];
      end
    
    FIELDS{1,ibox} = LOCAL_skel_proj(K_HSS{BOX.T_UP,ibox},K_HSS{BOX.J_UP,ibox},qloc);
  end
end

% Construct all incoming potentials:
FIELDS{2,2} = K_HSS{BOX.M_SIB,2}*FIELDS{1,3};
FIELDS{2,3} = K_HSS{BOX.M_SIB,3}*FIELDS{1,2};

for ibox = 2:nboxes
  if ( (K_HSS{BOX.C1,ibox}>0) || (K_HSS{BOX.C2,ibox}>0) ) % ibox has at least one son
      if (K_HSS{BOX.C1,ibox}>0)
          ison1   = K_HSS{BOX.C1,ibox};
          if (K_HSS{BOX.C2,ibox}>0)
              ison2   = K_HSS{BOX.C2,ibox};
              n1      = length(K_HSS{BOX.I_SKDN,ison1});
              n2      = length(K_HSS{BOX.I_SKDN,ison2});
              u_long  = LOCAL_skel_eval(K_HSS{BOX.T_DN,ibox},K_HSS{BOX.J_DN,ibox},FIELDS{2,ibox});
              FIELDS{2,ison1} = u_long(1:n1,:)        + K_HSS{BOX.M_SIB,ison1}*FIELDS{1,ison2};
              FIELDS{2,ison2} = u_long(n1 + (1:n2),:) + K_HSS{BOX.M_SIB,ison2}*FIELDS{1,ison1};
          else
              FIELDS{2,ison1} = LOCAL_skel_eval(K_HSS{BOX.T_DN,ibox},K_HSS{BOX.J_DN,ibox},FIELDS{2,ibox});
          end
      else
          ison2   = K_HSS{BOX.C2,ibox};
          FIELDS{2,ison2} = LOCAL_skel_eval(K_HSS{BOX.T_DN,ibox},K_HSS{BOX.J_DN,ibox},FIELDS{2,ibox});
      end
  end
end

% Construct the potential on all the leaves:
uu = zeros(size(qq));
for ibox = nboxes:(-1):2
  if ( (K_HSS{BOX.C1,ibox}<=0) && (K_HSS{BOX.C2,ibox}<=0) )
    ind       = dim*K_HSS{BOX.END1,ibox} - dim + (1:dim*K_HSS{BOX.END2,ibox});
    uu(ind,:) = LOCAL_skel_eval(K_HSS{BOX.T_DN,ibox},K_HSS{BOX.J_DN,ibox},FIELDS{2,ibox}) + ...
                K_HSS{BOX.M_SELF,ibox}*qq(ind,:);
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function that applies R = [I / T]
function qq_proj = LOCAL_skel_proj(T,J,qq)

if (size(T,2) == 0)
  qq_proj = qq(J,:);
else
  k       = size(T,1);
  qq_proj = qq(J(1:k),:) + T*qq(J((k+1):end),:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function that applies L = [I // T']
function uu_eval = LOCAL_skel_eval(T,J,uu)

uu_eval = zeros(length(J),size(uu,2));
if (size(T,2) == 0)
  uu_eval(J,:) = uu;
else
  k                       = size(T,1);
  uu_eval(J(1:k),:)       = uu;
  uu_eval(J((k+1):end),:) = T.'*uu;
end
