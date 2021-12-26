
//-------------------------------Test Richardson-----------------------------

//initialisation de la matrice A

function [A]=A_poisson1D(m)
  A = zeros(m, m);
  k=1
  while k <=m
    A(k, k) = 2;
    k=k+1
  end
  k=1
  while k <m
    A(k + 1, k) = -1;
    A(k, k + 1) = -1;
    k=k+1;
  end
endfunction

//initialisation de b

function [b] =b_poisson1D(m)
    b = zeros (m, 1);
    b(1) = -5;
    b(m) = 5;
endfunction

//------------------------Fonction test ------------------------
/*
cette fonction a pour but de test la methode de Richardson 

        itmax: nombre maximum de l'itération
        tol: tolérance de l'erreur
*/
//-------------initialisation de parametres d'entrées --------------------
m=3
itmax=100
tol=1e-8
x0=[1,1,1]'
alpha=0.25
    
function [x,relres,vec,k]=test_Richardson(m, itmax, tol,x0,alpha)

    A = A_poisson1D(m);
    disp(A)
    b =b_poisson1D(m)
    disp(b)
    [x,relres,vec,k]=Richardson_P1D(A,b,tol,itmax,x0,alpha)


endfunction
