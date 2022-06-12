n=input('What is the dimension n? ');
h=1/n;n=n-1;A=zeros(n);Q=eye(n);R=zeros(n);
for i=1:n
    A(i,i)=2/h^2;
end
for i=1:n-1
    A(i,i+1)=-1/h^2;
    A(i+1,i)=-1/h^2;
end
saveA=A;
for j=1:n
    q=A(:,j);
    for i=1:j-1
        R(i,j)=Q(:,i)'*A(:,j);
        q=q-R(i,j)*Q(:,i);
    end
    R(j,j)=norm(q);
    Q(:,j)=q/R(j,j);
end
for j=1:n
    for i=1:n
        if(abs(A(i,j)-Q(i,:)*R(:,j))>5*10^-13)
        fprintf('The computed solution seems to be wrong at %f \n',j);
        end
    end
end
disp('The Gram_Schmidt decomposition is as follows: ')
disp('Q')
disp(Q)
disp('R')
disp(R)
