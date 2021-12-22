//----------------------gauss_Seidel------------------------------

   /*Valeur retour:
        x: est la résolution du système
        k: est nombre d'intération à converger
        errl: l’erreur relative
    */
function [x, errl, k] = Gss_Seidel(A, b, tol, itmax)
    /*
        x0: valeur initiale
        itmax: nombre maximum de l'itération
        tol: tolérance de l'erreur
    */    
    n = size(A,"c");
    D = zeros(n, n);
    for i = 1:n
      D(i, i) = A(i, i);
    end
    E = zeros(n, n);
    for i = 2:n
      E(i, i - 1) = - A(i, i - 1);
	end
	DE = inv(D - E);
	x0 = zeros(n, 1);

	errl=zeros(itmax,1);
   
    res = b-A*x0;
    relres = norm(res)/norm(b);
    
    for k = 1:itmax-1 
		x = x0 + DE * res;
        x0=x
		r = (b - A * x0);

		res = b-A*x0;
        relres = norm(res)/norm(b);
        if relres<tol
           break
        end
        errl(k,1) = relres;
     end

endfunction
