format("e", 16);

//----------------initialisation des parametres d'entrées----------------
n = 40;
e = 10^(-10);
max_count = 10^6;

y_jacobi = zeros(n);
y_GSeidel = zeros(n);
y_jacobi_relres = zeros(n);
y_GSeidel_relres = zeros(n);

//-------------initialisation de la matrice A -------------------------

function [A] = A_poisson1D(n)

  A = zeros(n, n);
  k=1
  while k <=n
    A(k, k) = 2;
    k=k+1
  end
  k=1
  while k <n
    A(k + 1, k) = -1;
    A(k, k + 1) = -1;
    k=k+1;
  end
endfunction

//----------------------initialisation du vecteur b-----------------------

function [b] = b_poisson1D(n)
    
  b = zeros (n, 1);
  b(1) = -5;
  b(n) = 5;
endfunction

/*
Remarque : ici j'ai ajouté le parametre de sortie relres qui signifier l'erreur relative pour que je puisse faire les graphes et notre analyse pour la suite
*/

function [x, k, relres,errl] = jacobi(A,b,tol,itmax)
    /*
        x0: valeur initiale
        itmax: nombre maximum de l'itération
        tol: tolérance de l'erreur
    */    
    n = size(A,"c");
    x0 = zeros(n,1);
    errl=zeros(itmax,1);
   
    r = b-A*x0;
    relres = norm(r)/norm(b);
    D = (1.)./diag(A);
    
    for k = 1:itmax-1    
        x = x0 + D.*r;
        x0 = x;
        r = b-A*x0;
        relres = norm(r)/norm(b);
        if relres<tol
           break
        end
        errl(k,1) = relres;
    end    
endfunction

function [x, k, relres,errl] = gauss_seidel(A, b, tol, itmax)
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
	SB = inv(D - E);
	x0 = zeros(n, 1);

	errl=zeros(itmax,1);
   
    r = b-A*x0;
    relres = norm(r)/norm(b);
    
    for k = 1:itmax-1 
		x = x0 + SB * r;
        x0=x
		r = (b - A * x0);
        relres = norm(r)/norm(b);
        if relres<tol
           break
        end
        errl(k,1) = relres;
     end

endfunction

k=1
while k<=n

  A = A_poisson1D(k);
  b = b_poisson1D(k);

  [x_j, k_j, relres_j, errl_j] = jacobi(A, b, e, max_count);
  [x_g, k_g, relres_g, errl_g] = gauss_seidel(A, b, e, max_count);

  y_jacobi(k) = k_j;
  y_GSeidel(k) = k_g;

  y_jacobi_relres(k) = relres_j;
  y_GSeidel_relres(k) = relres_g;

  k=k+1
end


//-------------------------Taille matrice vs erreur relative---------------------

xtitle("Evolution de l erreur relative par rapport au taille de la matrice", "Taille de la matrice", "erreur relative");
plot([1:n], [log(y_jacobi_relres) log(y_GSeidel_relres)]);
legend(["jacobi" "gauss seidel"], 2);
xs2png(0, "img/number_of_iteration_relres.png");
clf();


//----------------------------------Erreur---------------------------------------

n = 10;
y_jacobi_e = zeros(n);
y_GSeidel_e = zeros(n);
y_jacobi_e_relres = zeros(n);
y_GSeidel_e_relres = zeros(n);

error_vec = [10^(-1), 10^(-2), 10^(-3), 10^(-4), 10^(-5), 10^(-6), 10^(-7), 10^(-8), 10^(-9), 10^(-10), 10^(-11), 10^(-12), 10^(-13), 10^(-14), 10^(-15)];

k=1;
while k<=12

  A = A_poisson1D(k);
  b = b_poisson1D(k);

  e = error_vec(k);

  [x_j, k_j, relres_j, errl_j] = jacobi(A, b, e, max_count);
  [x_g, k_g, relres_g, errl_g] = gauss_seidel(A, b, e, max_count);

  y_jacobi_e(k) = k_j;
  y_GSeidel_e(k) = k_g;

  y_jacobi_e_relres(k) = relres_j;
  y_GSeidel_e_relres(k) = relres_g;

  k=k+1
end


//------------------------Erreur relative vs Taille de la matrice ----------------

xtitle("L evolution de l ereeur relative", "Taille de la matrice", "Erreur relative");
plot([1:12], [log(y_jacobi_e_relres) log(y_GSeidel_e_relres)]);
legend(["Jacobi" "Gauss Seidel"], 2);
xs2png(0, "img/Erreur_relatives.png");
clf();

//-------------------------------------- Convergence-----------------------------
n = 10;
e = 10^(-8);
max_count = 600;

A = A_poisson1D(n);
b = b_poisson1D(n);
  
[x_j, k_j, relres_j, errl_j] = jacobi(A, b, e, max_count);
[x_g, k_g, relres_g, errl_g] = gauss_seidel(A, b, e, max_count);
disp(max(size(relres_j)));

xtitle("Convergence", "nombre d iteration", "Erreur");
plot(log10(errl_j), "g")
plot(log10(errl_g))
legend(["jacobi" "gauss seidel"], 1);
xs2png(0, "img/convergence.png");
clf();
