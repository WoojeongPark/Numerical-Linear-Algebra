function G = Cholesky_decompose(A)
% Assume A is symmetric, positive - definite
n = size(A,1);
G = eye(n);

% Cholesky decomposition
for j = 1 : n
    G(j,j) = sqrt(A(j,j) - dot(G(j, 1:j-1), G(j, 1:j-1)));
    for k = (j+1) : n
       G(k,j) = (A(k,j) - dot(G(k, 1:j-1), G(j, 1:j-1)))/G(j,j); 
    end
end

A_new = G*G';

% check solution
tolerance = 1e-10;
count = 0;
for i = 1 : n
   for j = 1 : n
       if(abs(A(i,j) - A_new(i,j)) > tolerance) 
           count = count + 1;
       end
   end
end

if(count == 0) 
    disp("correct decomposition");
else
    disp("incorrect decomposition");
end
end