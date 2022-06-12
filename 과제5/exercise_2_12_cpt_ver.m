h=input('how many meshes do you want? Please enter the # of partition h(=1/N): ');
n=1/h;n=n-1;
A=zeros(n,n);b=zeros(n,1);al=zeros(n,1);au=zeros(n,1);ad=zeros(n,1);
saveb=zeros(n,1);saveA=zeros(n,n);
for j=1:n %필요한 행렬 생성. 이 때, u(0)=0, u(1)=0을 가정하고, 후에 1을 더하여 해를 구한다.
    A(j,j)=2;
end
for j=1:n-1
    A(j,j+1)=-1;
    A(j+1,j)=-1;
end
for k=1:n
    b(k)=h^2*exp(sin(k/(n+1)));
end
saveA=A;
saveb=b;
for j=2:n %Thomas algorithm setting
    al(j)=-1;
    au(j-1)=-1;
    ad(j)=2;
end
ad(1)=2; 
for j=2:n %Thomas algorithm
    m=al(j)/ad(j-1); al(j)=m;
    ad(j)=ad(j)-m*au(j-1);
    b(j)=b(j)-m*b(j-1);
end
b(n)=b(n)/ad(n);%back substitution, overwriting on b
for j=n-1:-1:1 
    b(j)=(b(j)-au(j)*b(j+1))/ad(j);
end
for j=1:n %check the solution is correct;
    if(abs(saveb(j)-dot(saveA(j,:),b(:)))>5*10^-13)
    fprintf('The computed solution seems to be wrong at %f \n',j);
    end
end
disp('The solution x is as follows: ')
b=b+1; %최종적으로 boundary condition을 만족시키기 위해 1을 더한다. 초기값 설정
b(1)=1;b(n)=1;
disp(b)

range = linspace(0,1,n);
title('Exercise 2_28','fontsize',12,'fontname','arial');
xlabel('x');
ylabel('solution');
plot(range,b,'-o');
legend('n=128')