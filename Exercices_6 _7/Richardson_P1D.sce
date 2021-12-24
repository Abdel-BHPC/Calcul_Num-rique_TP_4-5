//----------------------Richardson------------------------------

   /*Valeur retour:
        x: est la résolution du système
        k: est nombre d'intération à converger
        relres: l’erreur relative
        vec: vecteur de l'erreur relative
    */
function[x,relres,vec,it]=Richardson_P1D(A,b,tol,itmax,x0,alpha)
    /*
        x0: valeur initiale
        itmax: nombre maximum de l'itération
        tol: tolérance de l'erreur
    */ 
    n=size(A,"r");
    vec = zeros(itmax,1);
    r=b-A*x0
    relres=norm(r)/norm(b);
    for k = 1:itmax-1 
        k=k+1;
        x=x0+alpha*r;
        x0=x;
        r=b-A*x0;
        relres=norm(r)/norm(b);
        if relres<tol
           break
        end
        vec(k)=relres;
    end
endfunction   

