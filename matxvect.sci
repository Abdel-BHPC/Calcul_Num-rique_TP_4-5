function [y] = matxvect(AA, IA, JA, x)
    n = size(JA,2);
    y = zeros(n-1, 1);
    for i=1:n-1
        for j=JA(i):JA(i+1)-1 
            y(i) = y(i) + AA(j)*x(IA(j))
        end
    end
endfunction

