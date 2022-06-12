n=input('What is the dimension n? ');
h=1/n;n=n-1;
A=zeros(n,n);b=zeros(n,1);c=zeros(n+1,1);x=zeros(n,1); %c must have n+1 dimension.
saveb=zeros(n,1);saveA=zeros(n,n);savec=zeros(n+1,1);
for k=1:n
    b(k)=h^2*(-(pi)*cos(2*pi*h*k)+(2*pi)*sin(pi*h*k));
    c(k)=sin(pi*h*(2*k-1)/2)+2;
end
c(n+1)=sin(pi*h*(2*n+1)/2)+2; %We need one more dimension.
for j=1:n %Matrix for solving the equation.
    A(j,j)=c(j+1)+c(j);
end
for j=1:n-1
    A(j,j+1)=-c(j+1);
    A(j+1,j)=-c(j+1);
end
saveA=A;saveb=b;savec=c;
for j=1:n-1 %Gauss elimination
    for k=j+1:n
        m=A(k,j)/A(j,j);
        A(k,j:n)=A(k,j:n)-m*A(j,j:n);
        b(k)=b(k)-m*b(j);
    end
end
for j=n:-1:1 %find solution(overwriting b) by back substitution.
    b(j)=(b(j)-dot(A(j,j+1:n),b(j+1:n)))/A(j,j);
end
for j=1:n
if(abs(saveb(j)-dot(saveA(j,1:n),b(1:n)))>5*10^-13)
    fprintf('The computed solution seems to be wrong at %f \n',j);
end
end
disp('The solution x is as follows: ')
disp(b)
range = linspace(0,1,n);
title('Exercise 2_23','fontsize',12,'fontname','arial');
xlabel('x');
ylabel('solution');
plot(range,b,'-o');
legend('n=100')

