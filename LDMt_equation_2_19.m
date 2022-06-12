n=input('What is the dimension n? ');
%WANT : solve -u''(x)=1
%When considering the initial value, u(0)=0, u(1)=1, we may start to solve
%-y''(x)=1 when y(x)=u(x)-x. Then initial values of y is : y(0)=y(1)=0; 
h=1/n;n=n-1;y=zeros(n,1);b=ones(n,1);
for i=1:n
    b(i,1)=h^2*b(i,1);
end
%Decompose tridiagonal matrix to LDM^t
A=zeros(n,n);D=zeros(n,n);L=zeros(n);M=zeros(n);
r=zeros(1,n);w=zeros(1,n);
for j=1:n
    A(j,j)=2;
end
for j=1:n-1
    A(j,j+1)=-1;
    A(j+1,j)=-1;
end
saveA=A;
for j=1:n
    for p=1:j-1
        r(p)=D(p,p)*M(j,p);
        w(p)=L(j,p)*D(p,p);
    end
    D(j,j)=A(j,j)-dot(w(1:j-1),M(j,1:j-1));
    for k=j+1:n
        L(k,j)=(A(k,j)-dot(L(k,1:j-1),r(1:j-1)))/D(j,j);
        M(k,j)=(A(k,j)-dot(w(1:j-1),M(k,1:j-1)))/D(j,j);
    end
end
for i=1:n
    L(i,i)=1;
    M(i,i)=1;
end
%Now, let's find the solution. Put (D M^t y)=y1, then the equation is equivalent to Ly1=b. Find y1 by forward substitution.
y1=zeros(n,1);
y1(1)= b(1)/L(1,1);
for j = 2 : n
y1(j) = (b(j) - dot(L(j,1:j-1), y1(1:j-1,1)))/L(j,j);
end
%Next, put (M^t y)=y2, then the equation is: Dy2=y1. Find y2. It is very
%easy since D is diagonal matrix.
y2=zeros(n,1);
for j=1:n
    y2(j)=y1(j)/D(j,j);
end
%At last, the equation is M^t * y=y2, and M^t is unit upper triangular matrix, so we must find y using back
%substitution.
M=M';
y(n) = y2(n); %M(n,n)=1, so we don't need to divide it.
for j = (n-1) : -1 : 1
y(j) = (y2(j) - dot(M(j,j+1:n), y(j+1:n,1)))/M(j,j);
end
%The final solution is u(x)=y(x)+x, and each y(i) in matlab means y(i/n)
%where n is the dimension. For n=n-1 above, we should code like;
u=zeros(n,1);
for i=1:n
    u(i)=y(i)+i/(n+1);
end
u
%Note that if n is sufficiently large, u(0) is almost 0 and u(1) is almost
%1.
range = linspace(0,1,n);
title('Exercise 2_19','fontsize',12,'fontname','arial');
xlabel('x');
ylabel('solution');
plot(range,u,'-o');
legend('n=64')