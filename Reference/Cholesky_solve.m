function x = Cholesky_solve(n)
%Initialization
A = zeros(n);
for j = 1 : n
    A(j,j) = 2; 
end

for j = 1 : (n-1)
   A(j,j+1) = -1; 
   A(j+1,j) = -1;
end
    
b = zeros(n,1);
b(1) = 1;

G = Cholesky_decompose(A);
% now we solve GG'x = b, G is lower triangular matrix
% first, forward substitution (Gy = b)
y = zeros(n,1);
y(1) = b(1)/G(1,1);
for j = 2 : n
    y(j) = (b(j) - dot(G(j,1:j-1), y(1:j-1)))/G(j,j);
end

% second, back substitution (Ux = b, U = G')
U = G';
x = zeros(n,1);
x(n) = y(n)/U(n,n);
for j = (n-1) : -1 : 1
     x(j) = (y(j) - dot(U(j,j+1:n), x(j+1:n,1)))/U(j,j);
end

% check solution
b_new = A*x;
tolerance = 1e-10;
count = 0;
for i = 1 : n
    if(abs(b(i) - b_new(i)) > tolerance) 
        count = count + 1;
    end
end

if(count == 0) 
    disp("correct solution");
else
    disp("incorrect solution");
end

end