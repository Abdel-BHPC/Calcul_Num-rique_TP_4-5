//----------------------Jacobi------------------------------

   /*Valeur retour:
        x: est la résolution du système
        k: est nombre d'intération à converger
        errl: l’erreur relative
    */
function [x, errl, k] = Jaccobi(A,b,tol,itmax)
    /*
        x0: valeur initiale
        itmax: nombre maximum de l'itération
        tol: tolérance de l'erreur
    */    
    n = size(A,"c");
    x0 = zeros(n,1);
    errl=zeros(itmax,1);
   
    res = b-A*x0;
    relres = norm(res)/norm(b);
    D = (1.)./diag(A);
    
    for k = 1:itmax-1    
        x = x0 + D.*res;
        x0 = x;
        res = b-A*x0;
        relres = norm(res)/norm(b);
        if relres<tol
           break
        end
        errl(k,1) = relres;
    end    
endfunction
