n=input('What is the dimension n? ');
h=1/n;n=n-1;A=zeros(n);Q=eye(n);R=zeros(n);b=zeros(n,1);x=zeros(n,1);
for i=1:n
    A(i,i)=2/h^2;
    b(i,1)=0;
end
b(1,1)=1;saveb=b;
for i=1:n-1
    A(i,i+1)=-1/h^2;
    A(i+1,i)=-1/h^2;
end
saveA=A;
for i=1:n
    h=house(A(i:n,i)); %householder function. 
    H=eye(n);
    H(i:n,i:n)=h;
    Q=Q*H;
    A=H*A;
end
R=A;
b=Q'*b;
for j=n:-1:1
    x(j,1)=(b(j,1)-dot(R(j,j+1:n),x(j+1:n,1)))/R(j,j);
end
for j=1:n
if(abs(saveb(j,1)-dot(saveA(j,1:n),x(1:n,1)))>5*10^-13)
    fprintf('The computed solution seems to be wrong at %f \n',j);
end
end
disp('The solution x is as follows: ')
disp(x)
