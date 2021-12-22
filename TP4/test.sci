
function[] = test()
A=[2,-1,0;-1,2,-1;0,-1,1] // matrice symetrique definie positive
[L,D,LT] = facldlt(A)
disp("L*D*LT =", L*D*LT)
disp("A =",A)
erreur = norm((A-L*D*LT) / A) 
disp("erreur = ",erreur)

endfunction
