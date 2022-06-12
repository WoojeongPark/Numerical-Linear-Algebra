A=zeros(7);b=zeros(7,1);
for k=1:10
    t(k)=k*pi/10;
    y(k)=cos(t(k));
end
for i=1:7
    for j=1:7
        ipower=14-i-j;
        for k=1:10
        A(i,j)=A(i,j)+t(k)^ipower;
        end
    end
end
for i=1:7
    for k=1:10
        b(i)=b(i)+t(k)^(7-i)*y(k);
    end
end
saveA=A;saveb=b;
G=zeros(7);count=0;
for i=1:7 %Cholesky decomposition
   G(i,i)=sqrt(A(i,i)-dot(G(i,1:i-1),G(i,1:i-1)));
   for j=(i+1):7
      G(j,i)=(A(j,i) -dot(G(j,1:i-1),G(i,1:i-1)))/G(i,i);
   end
end
for j=1:7 %Solve Gy=b where y=G'x
    b(j)=(b(j)-dot(G(j,1:j-1),b(1:j-1)))/G(j,j); %forward substitution
end
G=G';
for j=7:-1:1 %Solve G'x=y
    b(j)=(b(j)-dot(G(j,j+1:7),b(j+1:7)))/G(j,j); %backward substitution overwriting.
end
for j=1:7 %solution check
if(abs(saveb(j)-dot(saveA(j,1:7),b(1:7)))>5*10^-12)
    count=count+1;
end
end
if(count~=0)
    fprintf('The computed solution seems to be wrong at %f \n',j);
end
disp('The solution x is as follows: ')
disp(b)
range = linspace(0,pi,10);
title('Exercise 6_11','fontsize',12,'fontname','arial');
xlabel('t');
ylabel('solution');
plot(t,y,'-o');
hold on
x=linspace(0,pi,1000);poly=zeros(1000,1);
for i=1:1000
    poly(i,1)=b(1)*x(i)^6+b(2)*x(i)^5+b(3)*x(i)^4+b(4)*x(i)^3+b(5)*x(i)^2+b(6)*x(i)+b(7);
end
plot(x,poly);
hold off