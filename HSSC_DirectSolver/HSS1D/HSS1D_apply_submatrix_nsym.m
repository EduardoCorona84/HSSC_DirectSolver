function uu = HSS1D_apply_submatrix_nsym(K_HSS,Box_dn,Box_up,I_dn,I_up,qq,dim)
%
% This file is part of HSSC_DirectSolver
% Copyright (C) 2011-2013 Eduardo Corona, Per Gunnar Martinsson, Denis Zorin
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%    FUNCTION CALL:
%        uu = HSS1D_apply_submatrix_nsym(K_HSS,Box_dn,Box_up,I_dn,I_up,qq,dim)
%
%    DESCRIPTION:
%        This function applies a submatrix of a non-symmetric HSS1D matrix to a vector :
%            uu = K(I_dn,I_up) qq(I_up)
%
%    INPUT:
%        K denotes the matrix encoded by HSS1D structure K_HSS.
%        I_dn / I_up are boolean vectors that specify the row / column index sets of the submatrix
%            of K.
%        Box_dn / Box_up are the boolean index vectors specifying the boxes of K_HSS that
%            correspond to the submatrix. They can be precomputed as follows:
%            [B_dn,B_up] = HSS1D_apply_submatrix_boxes(K_HSS,I_dn,I_up);
%        qq is the vector to be left multiplied by the submatrix.
%        dim (int) is the dimension of the kernal.
%
%    OUTPUT:
%        uu is he vector resulting from the matrix-vector multiply
%            uu = K(I_dn,I_up) qq(I_up)
%


global BOX 

nboxes = size(K_HSS,2);
FIELDS = cell(2,nboxes);
q2 = size(qq,2);
I_srcup = cell(nboxes,1);
I_srcdn = I_srcup;

% Initialize the fields on the leaves
%  - Construct the outgoing fields.
%  - Set the incoming fields to zero.
for ibox = nboxes:(-1):2
    
  if ( (K_HSS{BOX.C1,ibox}<=0) && (K_HSS{BOX.C2,ibox}<=0) )
    ind = dim*K_HSS{BOX.END1,ibox} - 1 + (1:dim*K_HSS{BOX.END2,ibox});
    I_srcdn{ibox} = I_dn(ind);
    I_srcup{ibox} = I_up(ind);
    
    if Box_up(ibox) == true
        if (sum(I_srcup{ibox})==K_HSS{BOX.END2,ibox})
            FIELDS{1,ibox} = LOCAL_skel_proj(K_HSS{BOX.T_UP,ibox},K_HSS{BOX.J_UP,ibox},qq(ind,:));
        else
            FIELDS{1,ibox} = LOCAL_skel_proj_rest(K_HSS{BOX.T_UP,ibox},K_HSS{BOX.J_UP,ibox},I_srcup{ibox},qq(ind,:));
        end
    else
        FIELDS{1,ibox} = zeros(K_HSS{BOX.KSKEL,ibox},q2); 
    end
    FIELDS{2,ibox} = zeros(length(K_HSS{BOX.I_SKDN,ibox}),q2);
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
    
      if Box_up(ibox)==true
          FIELDS{1,ibox} = LOCAL_skel_proj(K_HSS{BOX.T_UP,ibox},K_HSS{BOX.J_UP,ibox},qloc);
      else
          FIELDS{1,ibox} = zeros(K_HSS{BOX.KSKEL,ibox},q2); 
      end
  end
end

% Construct all incoming potentials:
if Box_up(3)==true && Box_dn(2)==true
    if ( (K_HSS{BOX.C1,3}>0) || (K_HSS{BOX.C2,3}>0) ) || size(K_HSS{BOX.M_SIB,2},2)==K_HSS{BOX.KSKEL,3}
        FIELDS{2,2} = K_HSS{BOX.M_SIB,2}*FIELDS{1,3};
    else
        FIELDS{2,2} = K_HSS{BOX.M_SIB,2}(:,I_srcup{3})*FIELDS{1,3};
    end
else
    FIELDS{2,2} = zeros(length(K_HSS{BOX.I_SKUP,2}),q2);
end

if Box_up(2)==true && Box_dn(3)==true
    if ( (K_HSS{BOX.C1,2}>0) || (K_HSS{BOX.C2,2}>0) ) || size(K_HSS{BOX.M_SIB,3},2)==K_HSS{BOX.KSKEL,2}
        FIELDS{2,3} = K_HSS{BOX.M_SIB,3}*FIELDS{1,2};
    else
        FIELDS{2,3} = K_HSS{BOX.M_SIB,3}(:,I_srcup{2})*FIELDS{1,2};
        
    end
else
    FIELDS{2,3} = zeros(length(K_HSS{BOX.I_SKDN,3}),q2);
end

for ibox = 2:nboxes
  if ( (K_HSS{BOX.C1,ibox}>0) || (K_HSS{BOX.C2,ibox}>0) ) % ibox has at least one son
      if (K_HSS{BOX.C1,ibox}>0)
          ison1   = K_HSS{BOX.C1,ibox};
          if (K_HSS{BOX.C2,ibox}>0)
              ison2   = K_HSS{BOX.C2,ibox};
              n1      = length(K_HSS{BOX.I_SKDN,ison1});
              n2      = length(K_HSS{BOX.I_SKDN,ison2});
              
              if Box_dn(ibox)==true
                u_long  = LOCAL_skel_eval(K_HSS{BOX.T_DN,ibox},K_HSS{BOX.J_DN,ibox},FIELDS{2,ibox});
                if Box_up(ison1)==true && Box_up(ison2)==true
                    FIELDS{2,ison1} = u_long(1:n1,:)        + K_HSS{BOX.M_SIB,ison1}*FIELDS{1,ison2};
                    FIELDS{2,ison2} = u_long(n1 + (1:n2),:) + K_HSS{BOX.M_SIB,ison2}*FIELDS{1,ison1};
                elseif Box_up(ison2)==true
                    FIELDS{2,ison1} = u_long(1:n1,:)        + K_HSS{BOX.M_SIB,ison1}*FIELDS{1,ison2};
                    FIELDS{2,ison2} = u_long(n1 + (1:n2),:);
                elseif Box_up(ison1)==true
                    FIELDS{2,ison1} = u_long(1:n1,:);
                    FIELDS{2,ison2} = u_long(n1 + (1:n2),:) + K_HSS{BOX.M_SIB,ison2}*FIELDS{1,ison1};
                else
                    FIELDS{2,ison1} = u_long(1:n1,:);
                    FIELDS{2,ison2} = u_long(n1 + (1:n2),:);
                end
              else
                   FIELDS{2,ison1} = zeros(n1,q2); 
                   FIELDS{2,ison2} = zeros(n2,q2);
              end
          else
              if Box_dn(ison1)==true
                  FIELDS{2,ison1} = LOCAL_skel_eval(K_HSS{BOX.T_DN,ibox},K_HSS{BOX.J_DN,ibox},FIELDS{2,ibox});
              else
                  FIELDS{2,ison1} = zeros(length(K_HSS{BOX.I_SKDN,ison1}),q2); 
              end
          end
      else
          ison2   = K_HSS{BOX.C2,ibox};
          if Box_dn(ison2)==true
              FIELDS{2,ison2} = LOCAL_skel_eval(K_HSS{BOX.T_DN,ibox},K_HSS{BOX.J_DN,ibox},FIELDS{2,ibox});
          else
              FIELDS{2,ison2} = zeros(length(K_HSS{BOX.I_SKDN,ison2}),q2);  
          end
      end
  end
end

% Construct the potential on all the leaves:
uu = zeros(size(qq));

for ibox = nboxes:(-1):2
  if ( (K_HSS{BOX.C1,ibox}<=0) && (K_HSS{BOX.C2,ibox}<=0) )
    ind       = dim*K_HSS{BOX.END1,ibox} - dim + (1:dim*K_HSS{BOX.END2,ibox});
    
    if Box_dn(ibox)==true
        if (sum(I_srcdn{ibox})==K_HSS{BOX.END2,ibox})
            if Box_up(ibox)==true
            uu(ind,:) = LOCAL_skel_eval(K_HSS{BOX.T_DN,ibox},K_HSS{BOX.J_DN,ibox},FIELDS{2,ibox}) + ...
                    K_HSS{BOX.M_SELF,ibox}*qq(ind,:);
            else
                uu(ind,:) = LOCAL_skel_eval(K_HSS{BOX.T_DN,ibox},K_HSS{BOX.J_DN,ibox},FIELDS{2,ibox});
            end
        else
            if Box_up(ibox)==true
            uu(ind(I_srcdn{ibox}),:) = LOCAL_skel_eval_submatrix(K_HSS{BOX.T_DN,ibox},K_HSS{BOX.J_DN,ibox},I_srcdn{ibox},FIELDS{2,ibox}) + ...
                    K_HSS{BOX.M_SELF,ibox}(I_srcdn{ibox},I_srcup{ibox})*qq(ind(I_srcup{ibox}),:);                            
            else
                uu(ind(I_srcdn{ibox}),:) = LOCAL_skel_eval_submatrix(K_HSS{BOX.T_DN,ibox},K_HSS{BOX.J_DN,ibox},I_srcdn{ibox},FIELDS{2,ibox});
            end
        end
    else
        
    end
  end
end


uu = uu(I_dn,:); 

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function that applies R = [I / T]
function qq_proj = LOCAL_skel_proj(T,J,qq)

if (min(size(T)) == 0)
  qq_proj = qq(J,:);
else
  k       = size(T,1);
  qq_proj = qq(J(1:k),:) + T*qq(J((k+1):end),:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function that applies R = [I / T]
function qq_proj = LOCAL_skel_proj_rest(T,J,I,qq)

if (min(size(T)) == 0)
  qq_proj = qq(J,:);
  %qq_proj(~I,:) = 0;
else
  k       = size(T,1); R(:,J) = [eye(k) T];
  qq_proj = R(:,I)*qq(I,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function that applies L = [I // T']
function uu_eval = LOCAL_skel_eval(T,J,uu)

uu_eval = zeros(length(J),size(uu,2));
if (min(size(T)) == 0)
  uu_eval(J,:) = uu;
else
  k                       = size(T,1);
  uu_eval(J(1:k),:)       = uu;
  uu_eval(J((k+1):end),:) = T.'*uu;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function that applies L = [I // T']
function uu_eval = LOCAL_skel_eval_submatrix(T,J,I,uu)

uu_eval = zeros(length(J),size(uu,2));
if (min(size(T)) == 0)
  uu_eval(J,:) = uu;
  uu_eval = uu_eval(I,:);
else
  k = size(T,1); L(J,:) = [eye(k) ; T.']; 
  uu_eval = L(I,:)*uu;
end