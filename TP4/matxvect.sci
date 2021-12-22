function [y]= matxvect(AA, IA, JA, x)
    n = size(IA,2);
    y = zeros(n-1, 1);
    for i=1:n-1
        for j=IA(i):IA(i+1)-1 
            y(i) = y(i) + AA(j)*x(JA(j))
        end
    end
endfunction


