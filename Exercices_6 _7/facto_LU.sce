//-------------Factorisation LU d'une matrice tridiagonale------------
/*
Cette fonction realise la factorisation LU d'une matrice tridiagonale en se basant 
sur l'elemination de Gauss.

Containtes:
    A : est une matrice tridiagonale. 
    L : est une matrice triangulaire inférieure ayant des 1 sur la diagonale obtenue après la factorisation.
    U : est une matrice triangulaire supérieure obtenue après la factorisation.
    
Valeur de retour :
    L,U
*/

function [L,U] =LU_tridiagonale(A)

n=size(A,"r");
for k=1:n-1
   for i=k+1:n
     A(i,k)=A(i,k)/A(k,k);
     for j=k+1:n
        A(i,j)=A(i,j)-A(i,k)*A(k,j);
     end;
   end;
end;
U=triu(A);
L=A-U
for i = 1:n
     L(i,i)=1;
end

endfunction;
