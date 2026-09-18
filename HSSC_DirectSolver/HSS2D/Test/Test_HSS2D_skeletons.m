function Err = test_HSS2D_skeletons(TREE,params)

depth = TREE.depth; 
X = params.X_source; 
h = params.h; 
layers = params.layers; 
Err = zeros(1,depth+1);

if params.transinv == 0
Err_box = zeros(1,length(TREE.BOX)); 
end

for lev=1:depth 
    if params.transinv == 0
        numlev = TREE.numlev(lev+1); 
    else
        numlev = 1; 
    end
    
    for i=1:numlev
    BOX = TREE.box_numbers(lev+1,i);
    cB = TREE.BOX(BOX).cent; 
    I_src = TREE.BOX(BOX).I_src; 
    I_sk = TREE.BOX(BOX).I_sk; X_sk = X(I_sk,:);  
    
    if params.transinv == 1
        J = TREE.LEV(lev+1).J; 
        T = TREE.LEV(lev+1).T; 
        k = TREE.LEV(lev+1).k;  
    else
        J = TREE.BOX(BOX).J; 
        k = TREE.BOX(BOX).k; 
        if params.sym == 1
            T = TREE.BOX(BOX).T; 
        else
            T = TREE.BOX(BOX).T_up; 
        end
    end
    
    I_rs = I_src(J(k+1:end)); X_rs = X(I_rs,:); 
    
    if mod(lev,2) == 0
       scl = lev/2;
       colrad = 2^(-scl); 
       rowrad = colrad; 
    else
       scl = (lev-1)/2;
       rowrad = 2^(-scl); 
       colrad = rowrad/2;
    end
    
    [Xp,Xp2] = build_box_proxy(cB,colrad,rowrad,2*layers,h);
    Xp = [Xp ; Xp2];
    
    if size(X_rs,1)>0
        if params.transinv == 1
            if ~isempty(T.V)
                Err(lev+1) = norm(Kernel_Eval(Xp,X_rs,params) - Kernel_Eval(Xp,X_sk,params)*T.U*T.V.')/norm(Kernel_Eval(Xp,X_rs,params));
            else
                Err(lev+1) = norm(Kernel_Eval(Xp,X_rs,params) - Kernel_Eval(Xp,X_sk,params)*T.U)/norm(Kernel_Eval(Xp,X_rs,params));
            end
        else
            if ~isempty(T.V)
                Err_box(BOX) = norm(Kernel_Eval(Xp,X_rs,params) - Kernel_Eval(Xp,X_sk,params)*T.U*T.V.')/norm(Kernel_Eval(Xp,X_rs,params));
            else
                Err_box(BOX) = norm(Kernel_Eval(Xp,X_rs,params) - Kernel_Eval(Xp,X_sk,params)*T.U)/norm(Kernel_Eval(Xp,X_rs,params));
            end
            
            Err(lev+1) = Err(lev+1) + Err_box(BOX); 
        end
    end
    end
    
    Err(lev+1) = Err(lev+1)/numlev; 
end