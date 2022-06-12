n=input('What is the dimension n? ');
col=zeros(n,1);row=zeros(1,n);
perm=eye(n); %Permutation matrix exchanging rows.
A=zeros(n,n); AI=zeros(n,n);
check=zeros(n,n);saveI=eye(n);
for j=1:n
    A(j,j)=2;
end
for j=1:n-1
    A(j,j+1)=-1;
    A(j+1,j)=-1;
end
saveA=A;
for j=1:n
    [x,ksave]=max(abs(A(j:n,j))); %Partial pivoting
    k=ksave+j-1;
    perm([j,k],:)=perm([k,j],:); %Row exchange in the permutation matrix
    A([j,k],:)=A([k,j],:); 
    m=1/A(j,j); %Start Gauss-Jordran inverse matrix algorithm
    A(j,j)=m;
    col=m*A(:,j);
    row=A(j,:);
    A=A-col*row;
    A(:,j)=col;
    A(j,:)=-1*m*row;
    A(j,j)=m;
end
for i=1:n
    for j=1:n
        AI(i,j)=A(i,:)*perm(:,j);
    end
end %undo the order of their rows.
disp(A)
for i=1:n %check whether the inverse is correct within the tolerance 10e-10
    for j=1:n
        check(i,j)=saveA(i,:)*AI(:,j);
        if(abs(check(i,j)-saveI(i,j))>10^-10)
            disp("This inverse is incorrect")
        end
    end
end