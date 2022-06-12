n=input('What is the dimension n? ');
h=1/n;n=n-1;nn=n^2;
A=zeros(nn,nn);L=eye(nn,nn);
b=zeros(nn,1);u=zeros(nn,1); %v=u-1, then the initial values corresponds to 1.
for i=1:nn
    A(i,i)=4;
    b(i)=1*h^2;
end
for i=0:n-1
    for j=1:n-1
        A(n*i+j,n*i+j+1)=-1;
        A(n*i+j+1,n*i+j)=-1;
    end
end
for i=1:nn-n
    A(i,i+n)=-1;
    A(i+n,i)=-1;
end
saveA=A;
saveb=b;
for j=1:nn-1 %LU-decomposition without pivoting(overwriting A=U)
    ml=min(j+n,nn);
    L(j+1:ml,j)=A(j+1:ml,j)/A(j,j);
    for k=j+1:ml
        mu=min(j+n,nn);
        A(k,j:mu)=A(k,j:mu)-L(k,j)*A(j,j:mu);
    end
end
disp("The LU decomposition is: ")
!disp(L);
!disp(A);
y=zeros(nn,1);%forward elimination for solving Ly=b
y(1)=b(1);
for j=2:nn
    k=max(1,j-n);
    y(j)=b(j)-dot(L(j,k:j-1),y(k:j-1,1));
end
x=zeros(nn,1);%back substitution Ux=y
for j=nn:-1:1
    k=min(j+n,nn);
    x(j)=(y(j)-dot(A(j,j+1:k),x(j+1:k,1)))/A(j,j);
end
for j=1:nn %check the solution is correct;
if(abs(saveb(j)-dot(saveA(j,:),x(:)))>5*10^-13)
    fprintf('The computed solution seems to be wrong at %f \n',j);
end
end
x=x+1; %At last, for finding u=x+1, where v is the solution we find under the assumption. So, we add 1 !
n=n+2; %When we divide with meshes n, we get n+1 points. Don't forget that we started this program with n=n-1.
u=zeros(n); 
u(1,:)=1;u(n,:)=1;u(:,1)=1;u(:,n)=1; %boundary condition
for i=2:n-1
    for j=2:n-1
       u(i,j)=x((n-2)*(i-2)+j-1); %rearrange the solution from vector x to matrix u(i,j)
    end   
end
%Let's check plot!
[X,Y]=meshgrid(0:h:1);
surf(X,Y,u);
xlabel('x');
ylabel('y');
zlabel('u(x,y)');
title('2.25(n=100)');
    