n=input('What is the dimension n? ');
A=zeros(n,n);
matx=eye(n); %x1, x2, ..., xn이 column 단위로 저장된 행렬. 즉 xi=matx(:,i)
check=zeros(n,n);
saveI=matx;
for j=1:n
    A(j,j)=2;
end
for j=1:n-1
    A(j,j+1)=-1;
    A(j+1,j)=-1;
end
saveA=A;
for j=1:n-1
    for k=j+1:n
        m=A(k,j)/A(j,j);
        A(k,j:n)=A(k,j:n)-m*A(j,j:n);
        matx(k,:)=matx(k,:)-m*matx(j,:); %Gauss elimination
    end
end
for j=1:n
for i=n:-1:1
    matx(i,j)=(matx(i,j)-dot(A(i,i+1:n),matx(i+1:n,j)))/A(i,i);
end %find xj=matx(:,j) for j=1:n
end
disp('The solution A^-1 is as follows: ')
disp(matx)
for k=1:n
    for j=1:n
        check(i,j)=saveA(i,:)*matx(:,j);
    end
end
if(abs(check(i,j)-saveI(i,j)>10^-10))
    disp("This inverse is incorrect")
else
    disp("This inverse is correct almost everywhere")
end
    