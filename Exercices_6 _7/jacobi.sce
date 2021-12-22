//----------------------Jacobi------------------------------

/*Valeur retour:
    x: est la résolution du système
    k: est nombre d'intération à converger
    errl: l’erreur relative
*/
function [x,errl,k]=jacobi(A,b,itmax,tol,x0)
    
    /*
        x0: valeur initiale
        Kmax: nombre maximum de l'itération
        tol: tolérance de l'erreur
    */    
   n=size(A,"c");
   x=zeros(n,1);
   errl=zeros(itmax,1);
   
   res = b-A*x0;
   relres = norm(res)/norm(b);
   
   for k = 1:itmax-1
       x(1,1)=(1/A(1,1))*(b(1)-A(1,2)*x0(2));
       x(n,1)=(1/A(n,n))*(b(n)-A(n,n-1)*x0(n-1));
       for i=2:n-1
           x(i,1)=(1/A(i,i))*(b(i)-A(i,i-1)*x0(i-1)-A(i,i+1)*x0(i+1))
       end
       res = b-A*x0;
       relres = norm(res)/norm(b);
       if relres<tol
           break
       end
       errl(k,1)=relres;
       x0=x    
   end
endfunction
