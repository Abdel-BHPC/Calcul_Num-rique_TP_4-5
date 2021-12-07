function [L,D,LT]=facldlt(A)
  
  n=size(A,1);
  
  LT=zeros(n,n);
  L=zeros(n,n);
  //L = tril(A)
  for i = 1:n
    L(i,i)=1;
  end
  D = zeros(n,1);
  //v = zeros(n,1);
  for i=1:n
     D(i)=A(i,i)-((L(i,1:i)^2)*D(1:i));
     for j=i+1:n
        
        L(j,i)=(A(j,i)-L(j,1:i)*L(i,1:i)*D(1:i))/D(i);
     end
  end


  //D=diag(d)
  LT = L'
  //for j=1:n
   //  for i=1:j-1
    //   v(i)=L(j:i)*D(i)
    // end         
    // D(j)=A(j,j)-(L(j,1:j-1)*v(1:j-1));
    
  //   L(j+1:n,j) =(A(j+1:n,j)-(L(j+1:n,1:j-1)*v(1:j-1))/D(j)
    // end
 // end

endfunction
