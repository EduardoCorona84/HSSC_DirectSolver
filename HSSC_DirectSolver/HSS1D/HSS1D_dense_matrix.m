function Kd = HSS1D_dense_matrix(K_HSS,sym)
%
% This file is part of HSSDirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%     This code outputs the full, dense matrix Kd encoded by K_HSS in HSS1D form.
%


global BOX 

nboxes = size(K_HSS,2);
N = K_HSS{7,1};
Kd = zeros(N,N); 

% Construct the potential on all the leaves:
for ibox = nboxes:(-1):1  
    ison1 = K_HSS{BOX.C1,ibox}; ison2 = K_HSS{BOX.C2,ibox};
    % Leaf boxes: 
  if ( ison1<=0 && ison2<=0 )
    ind       = K_HSS{BOX.END1,ibox} - 1 + (1:K_HSS{BOX.END2,ibox});
    % Self interaction matrices
    Kd(ind,ind) = K_HSS{BOX.M_SELF,ibox};  
    
    % Store interpolation operators
    Rhat = zeros(K_HSS{BOX.KSKEL},K_HSS{BOX.NSKEL}); 
    Rhat(:,K_HSS{BOX.J_UP,ibox}) = [eye(K_HSS{BOX.KSKEL,ibox}) K_HSS{BOX.T_UP,ibox}]; 
    
    if sym
        Lhat = Rhat.'; 
    else
        Lhat = zeros(size(Rhat.')); 
        Lhat(K_HSS{BOX.J_DN,ibox},:) = [eye(K_HSS{BOX.KSKEL,ibox}) ; K_HSS{BOX.T_DN,ibox}.'];
    end
    
  elseif ( ison1>0 && ison2>0 )
      ind1 = K_HSS{BOX.END1,ison1} - 1 + (1:K_HSS{BOX.END2,ison1});
      ind2 = K_HSS{BOX.END1,ison2} - 1 + (1:K_HSS{BOX.END2,ison2});
      
      % Store sibling interactions
        Kd(ind1,ind2) = K_HSS{BOX.LAMBDA,ison1}*K_HSS{BOX.M_SIB,ison1}*K_HSS{BOX.RHO,ison2}; 
      if sym
        Kd(ind2,ind1) = Kd(ind1,ind2).';      
      else
        Kd(ind2,ind1) = K_HSS{BOX.LAMBDA,ison2}*K_HSS{BOX.M_SIB,ison2}*K_HSS{BOX.RHO,ison1}; 
      end
      
      n1 = size(K_HSS{BOX.RHO,ison1},2); n2 = size(K_HSS{BOX.RHO,ison2},2); 
      R_child = [K_HSS{BOX.RHO,ison1} zeros(K_HSS{BOX.KSKEL,ison1},n2) ; ...
                 zeros(K_HSS{BOX.KSKEL,ison2},n1)  K_HSS{BOX.RHO,ison2}]; 
      Rhat    =  LOCAL_skel_proj(K_HSS{BOX.T_UP,ibox},K_HSS{BOX.J_UP,ibox},R_child);  
      
      if sym
          Lhat = Rhat.'; 
      else
          n1 = size(K_HSS{BOX.LAMBDA,ison1},1); n2 = size(K_HSS{BOX.LAMBDA,ison2},1);
          
          L_child = [K_HSS{BOX.LAMBDA,ison1} zeros(n2,K_HSS{BOX.KSKEL,ison2}) ; ...
                 zeros(n1,K_HSS{BOX.KSKEL,ison1})  K_HSS{BOX.LAMBDA,ison2}]; 
          Lhat    = LOCAL_skel_proj(K_HSS{BOX.T_UP,ibox},K_HSS{BOX.J_UP,ibox},L_child.').'; 
          %LOCAL_skel_eval(K_HSS{BOX.T_DN,ibox},K_HSS{BOX.J_DN,ibox},L_child);  
      end
  else
      ison = max(ison1,ison2);  
      Rhat    =  LOCAL_skel_proj(K_HSS{BOX.T_UP,ibox},K_HSS{BOX.J_UP,ibox},K_HSS{BOX.RHO,ison});  
      
      if sym
          Lhat = Rhat.'; 
      else
          Lhat    =  LOCAL_skel_eval(K_HSS{BOX.T_DN,ibox},K_HSS{BOX.J_DN,ibox},K_HSS{BOX.LAMBDA,ison2});    
      end
      
  end
  
  K_HSS{BOX.RHO,ibox}    = Rhat; 
  K_HSS{BOX.LAMBDA,ibox} = Lhat; 
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
  size(T)
  size(uu)
  k                       = size(T,1);
  uu_eval(J(1:k),:)       = uu;
  uu_eval(J((k+1):end),:) = T.'*uu;
end
