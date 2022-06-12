n=input('What is the dimension n? ');
n=n-1;
A=zeros(n,n);b=zeros(n,1); G=zeros(n,n);
saveb=zeros(n,1);saveA=zeros(n,n); count=0;
for j=1:n
    A(j,j)=2;
end
for j=1:n-1
    A(j,j+1)=-1;
    A(j+1,j)=-1;
end
b(1)=1;
saveA=A;
saveb=b;
for i=1:n %Cholesky decomposition
   G(i,i)=sqrt(A(i,i)-dot(G(i,1:i-1),G(i,1:i-1)));
   for j=(i+1):n
      G(j,i)=(A(j,i) -dot(G(j,1:i-1),G(i,1:i-1)))/G(i,i);
   end
end
disp("The matrix G after decomposition is: ")
G
for j=1:n %Solve Gy=b where y=G'x
    b(j)=(b(j)-dot(G(j,1:j-1),b(1:j-1)))/G(j,j); %forward substitution
end
G=G';
for j=n:-1:1 %Solve G'x=y
    b(j)=(b(j)-dot(G(j,j+1:n),b(j+1:n)))/G(j,j); %backward substitution overwriting.
end
for j=1:n %solution check
if(abs(saveb(j)-dot(saveA(j,1:n),b(1:n)))>5*10^-13)
    count=count+1;
end
end
if(count~=0)
    fprintf('The computed solution seems to be wrong at %f \n',j);
end
disp('The solution x is as follows: ')
disp(b)
