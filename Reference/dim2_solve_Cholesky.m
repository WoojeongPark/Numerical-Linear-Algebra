function dim2_solve_Cholesky(N)
% we will solve 2-dim elliptic problem using finite difference method
% -(¥Äu(x,y)) = f(x,y) and B.C : u(x,y) = g(x,y)
% f(x,y) = 1, g(x,y) = 1
% In this method, we will use Cholesky decomposition method
n = N-1;
h = 1/N;

% we solve linear system Ax = b
% get A
A = zeros(n^2);
for i = 0 : (n-1)
   index = (n*i); 
   for j = 1 : n
       A(index + j, index + j) = 4;
   end
    
   for j = 1 : (n-1)
      A(index + j, index + (j+1)) = -1;
      A(index + (j+1), index + j) = -1;
   end
end

for j = 1 : (n^2 - n)
   A(j, j + n) = -1; 
end

for j = 1 : (n^2 - n)
   A(j + n, j) = -1; 
end

% now Cholesky decomposition
bw = n;
G = zeros(n^2);

for j = 1 : n^2
    mj = max(j-bw, 1);
    G(j,j) = sqrt(A(j,j) -  dot(G(j, mj:(j-1)), G(j, mj:(j-1))));
    
    for k = (j+1) : min(n^2, j+bw)
       mk = max(k-bw, 1);
       G(k,j) = (A(k,j) - dot(G(k, mk:(j-1)), G(j, mk:(j-1))))/G(j,j);
    end
end
    
% check decomposition
tolerance = 1e-10;
count = 0;
A_new = G*G';

for i = 1 : n^2
    for j = 1 : n^2
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
    
b = h^2 * ones(n^2);
% now we solve GG'x = b    
% first, Gy = b 
y = zeros(n^2, 1);
y(1) = b(1)/G(1,1);
for j = 2 : n^2
   k = max(1, j - bw); 
   y(j) = (b(j) - dot(G(j, k : (j-1)), y(k : (j-1))))/G(j,j);
end

% second, G'x = y, U = G'
U = G';
x = zeros(n^2, 1);
x(n^2) = y(n^2)/U(n^2,n^2);
for j = (n^2-1) : -1 : 1    
    k = min(j + bw, n^2);
    x(j) = (y(j) - dot(U(j, (j+1) : k), x((j+ 1) : k)))/U(j,j);
end

% check solution
count = 0;
b_new = A*x;
for i = 1 : n^2
   if(abs(b(i) - b_new(i) > tolerance)) 
       count = count + 1;
   end
end

if(count == 0) 
    disp("correct solution");
else
    disp("incorrect solution");
end    

% for boundary condition : g(x,y) = 1
x = (x + 1); 
% make solution to matrix form
Z = zeros(N+1);
Z(1,:) = ones(1,N+1);
Z(N+1,:) = ones(1,N+1);
for i = 2 : N
    Z(i,1) = 1;
    Z(i,N+1) = 1;
    for j = 2 : N
       Z(i,j) = x(n*(i-2) + (j-1)); 
    end   
end

% plot solution
[X,Y] = meshgrid(0 : h: 1);
surf(X, Y, Z);
xlabel('x');
ylabel('y');
zlabel('u(x,y)');
title('Exercise 2.27');
end