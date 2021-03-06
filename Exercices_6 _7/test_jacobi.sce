
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
cette fonction a pour but de comparer les deux methodes Gauss-Seidel et Jacobi en terme d'erreur

        itmax: nombre maximum de l'itération
        tol: tolérance de l'erreur
*/
    
function test_jacobi (m, itmax, tol)

    A = A_poisson1D(m);
  
    b =b_poisson1D(m)

    [x, errl1, k] = Jaccobi(A,b, tol, itmax);
    [x, errl2, k] = Gss_Seidel(A,b, tol, itmax);
    step = [1:itmax];
    xtitle ( "Résidu relative");
    plot ( step, errl1, "r-");
    plot ( step, errl2, "b-");
    legend ( 'Jacobi', 'Gauss-Seidel');

endfunction
