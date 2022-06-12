n=input('What is the dimension n? ');
h=1/n;n=n-1;nn=n^2;
A=zeros(nn,nn);L=eye(nn,nn);
b=zeros(nn,1);u=zeros(nn,1); %v=u-1, then the initial values corresponds to 1.
for i=1:nn
    A(i,i)=4;
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
for i=1:n
    for j=1:n
        b(n*(i-1)+j)=2*pi^2*sin(pi*j/(n+1))*sin(pi*i/(n+1));
    end
end
saveA=A;
saveb=b;
G=zeros(nn);
for j=1:nn %Cholesky decomposition
    mj=max(j-n,1);
    G(j,j)=sqrt(A(j,j)-dot(G(j,mj:j-1),G(j,mj:j-1)));
   for k=j+1:min(nn,j+n)
      mk=max(k-n,1);
      G(k,j)=(A(k,j)-dot(G(k,mk:j-1),G(j,mk:j-1)))/G(j,j);
   end
end
y=zeros(nn,1);%forward elimination for solving Gy=b
y(1)=b(1)/G(1,1);
for j=2:nn
    k=max(1,j-n);
    y(j)=(b(j)-dot(G(j,k:j-1),y(k:j-1,1)))/G(j,j);
end
x=zeros(nn,1);%back substitution G'x=y
H=G';
x(nn)=y(nn)/H(nn,nn);
for j=(nn-1):-1:1
    k=min(j+n,nn);
    x(j)=(y(j)-dot(H(j,j+1:k),x(j+1:k,1)))/H(j,j);
end
for j=1:nn %check the solution is correct;
if(abs(saveb(j)-dot(saveA(j,:),x(:)))>5*10^-10)
    fprintf('The computed solution seems to be wrong at %f \n',j);
end
end
n=n+2; %When we divide with meshes n, we get n+1 points. Don't forget that we started this program with n=n-1.
u=zeros(n); 
u(1,:)=0;u(n,:)=0;u(:,1)=0;u(:,n)=0; %boundary condition
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
title('2.31(h=1/100,Cholesky)');
    