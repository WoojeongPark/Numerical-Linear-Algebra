n=input('What is the dimension n? ');
h=1/n;n=n-1;nn=n^2;
A=zeros(nn,nn);L=eye(nn,nn);
b=zeros(nn,1);u=zeros(nn,1); 
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
for j=1:nn-1
    [x,jsave]=max(abs(A(j:min(j+n,nn),j)));
    js=jsave+j-1;
    A([j,js],j:min(j+2*n,nn))=A([js,j],j:min(j+2*n,nn));
    b([j,js])=b([js,j]);
    for k=j+1:min(j+n,nn)
        m=A(k,j)/A(j,j);
        A(k,j:min(j+2*n,nn))=A(k,j:min(j+2*n,nn))-m*A(j,j:min(j+2*n,nn));
        b(k)=b(k)-m*b(j);
    end
end
x=zeros(nn,1);%back substitution Ux=b' (b' : partial pivoting ÀÌÈÄÀÇ b)
x(nn)=b(nn)/A(nn,nn);
for j=(nn-1):-1:1
    k=min(j+n,nn);
    x(j)=(b(j)-dot(A(j,j+1:k),x(j+1:k,1)))/A(j,j);
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
title('2.31(h=1/100,LU pivoting)');