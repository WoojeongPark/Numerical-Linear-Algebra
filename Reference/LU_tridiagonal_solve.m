function LU_tridiagonal_solve(N)
% we solve -u''(x) = exp(sin(x)) in (0,1) with B.C : u(0) = u(1) = 1
% we will solve this system using linear system Ax = b

h = 1/N;
n = N-1;
% set A
A = zeros(n);
for j = 1 : n
    A(j,j) = 2;
end

for j = 1 : (n-1) 
    A(j,j+1) = -1;
    A(j+1,j) = -1;
end

% set b
b = zeros(n,1);
for j = 1 : n
   b(j) = h^2 * (exp(sin(h*j))); 
end

% LU decomposition of tridiagonal matrix A
U = eye(n);
L = eye(n);

for j = 1 : (n)
   if(j == 1)
       U(1,1) = A(1,1);
       L(2,1) = A(2,1)/U(1,1);
   end
    
   if(j > 1 && j < n)
       U(j-1,j) = A(j-1,j);
       U(j,j) = A(j,j) - L(j,j-1)*U(j-1,j);
       L(j+1,j) = A(j+1,j)/U(j,j);
   end
   
   if(j == n)
      U(n-1,n) = A(n-1,n); 
      U(n,n) = A(n,n) - L(n,n-1)*U(n-1,n);
   end
end

% now we solve LUx = b
% first Ly = b
y = zeros(n,1);
y(1) = b(1)/L(1,1);
for j = 2 : n
    y(j) = (b(j) - dot(L(j,1:j-1), y(1:j-1)))/L(j,j);
end

% second Ux = y
x = zeros(n,1);
x(n) = y(n)/U(n,n);
for j = (n-1) : -1 : 1
     x(j) = (y(j) - dot(U(j,j+1:n), x(j+1:n,1)))/U(j,j);
end

x = x + 1;
u = zeros(N+1);
u(1) = 1;
u(N+1) = 1;

for j = 2 : N
    u(j) = x(j-1);
end

X = linspace(0, 1, N+1);
plot(X, u, '-o');
ylim([1,1.25]);
xlabel('x');
ylabel('u(x)');
title('Ex2.28');
end