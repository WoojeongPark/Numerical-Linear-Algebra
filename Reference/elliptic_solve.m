function elliptic_solve(N)
n = N-1;
h = 1/N;

% we will use result of Exercise 2.22
% c(x) = sin(pi*x) + 2, f(x) = -(pi) * cos(2*pi*x) + (2*pi) * sin(pi*x)
% solve -(c(x)*u'(x))' = f(x), u(0) = u(1) = 0

% we get linear system Ax = b for the elliptic differential equation
% set A
A = eye(n);
for j = 1 : n
    xjl = j*h - (h/2);
    xjr = j*h + (h/2);
    A(j,j) = (sin(pi*xjl) + 2) + (sin(pi*xjr) + 2);  
end

for j = 1 : (n-1)
    xjr = (j*h) + (h/2);
    A(j,j+1) = -(sin(pi*xjr) + 2);
    A(j+1,j) = -(sin(pi*xjr) + 2);
end

b = zeros(n,1);
for j = 1 : n
    xj = (j*h);
    b(j) = h^2 * ( -(pi) * cos(2 * pi * xj) + (2*pi) * sin(pi * xj));
end

% now we can use Cholesky decomposition for this linear system
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

% get real u
u = zeros(N+1,1);
u(1) = 0;
u(N+1) = 0;
for j = 1 : n
   u(j+1) = x(j); 
end

% we know that real solution is (1/pi)*sin(pi*x)
x = linspace(0, 1, N+1);
plot(x, u, '-ro');
hold on;
plot(x, (1/pi)*sin(pi*x));
xlabel('x');
ylabel('u(x)');
title('Ex2.23');
hold off;
end