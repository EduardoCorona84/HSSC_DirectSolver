function I_bd = Find_Interface(X,cB,lev,h,lay,type)
% Produces a boolean array which tags boundary points as true 
% and interface / non-boundary points as false.
% X is a set of skeleton points on lay boundary layers inside box with
% center cB at level lev. 
% type determines where the interface with its sibling box lies. 

% Center points X
CB = repmat(cB,size(X,1),1);
XC = X - CB;

% Dimensions
X1max = max(XC(:,1)); 
X2max = max(XC(:,2));
Lsize = lay*h - h/2;

k = size(X,1); 

% Separate in North, South, East and West
if mod(lev,2) == 1
    IN = (XC(:,2)>=(X2max - Lsize));
    IS = (XC(:,2)<=-(X2max - Lsize));
    IE = (abs(XC(:,2))<(X2max - Lsize) & (XC(:,1)>=(X1max - Lsize)));
    IW = (abs(XC(:,2))<(X2max - Lsize) & (XC(:,1)<=-(X1max - Lsize)));
else
    IN = (XC(:,2)>=(X2max - Lsize)) & (abs(XC(:,1))<(X1max - Lsize));
    IS = (XC(:,2)<=-(X2max - Lsize)) & (abs(XC(:,1))<(X1max - Lsize));
    IE = (XC(:,1)>=(X1max - Lsize));
    IW = (XC(:,1)<=-(X1max - Lsize));
end

I_bd = true(1,k); 

if mod(lev,2) == 0
    if strcmp(type,'north') || strcmp(type,'west')
        I_bd(IS) = false; 
    elseif strcmp(type,'east') || strcmp(type,'south')
        I_bd(IN) = false; 
    end
else
    if strcmp(type,'north') || strcmp(type,'west')
        I_bd(IE) = false; 
    elseif strcmp(type,'south') || strcmp(type,'east')
        I_bd(IW) = false; 
    end
end