n=input('What is the dimension n? ');
A=zeros(n,n);D=zeros(n,n);L=zeros(n);M=zeros(n);
r=zeros(1,n);w=zeros(1,n);
check1=zeros(n,n);check2=zeros(n,n);
for j=1:n
    A(j,j)=2;
end
for j=1:n-1
    A(j,j+1)=-1;
    A(j+1,j)=-1;
end
saveA=A;
for j=1:n
    for p=1:j-1
        r(p)=D(p,p)*M(j,p);
        w(p)=L(j,p)*D(p,p);
    end
    D(j,j)=A(j,j)-dot(w(1:j-1),M(j,1:j-1));
    for k=j+1:n
        L(k,j)=(A(k,j)-dot(L(k,1:j-1),r(1:j-1)))/D(j,j);
        M(k,j)=(A(k,j)-dot(w(1:j-1),M(k,1:j-1)))/D(j,j);
    end
end
for i=1:n
    L(i,i)=1;
    M(i,i)=1;
end
disp("The LDM^t decomposition result is: ")
L=L
D=D
M=M
for i=1:n %check whether the inverse is correct within the tolerance 10e-10
    for j=1:n
        check1(i,j)=L(i,:)*D(:,j);
    end
end
M=M';
for i=1:n
    for j=1:n
        check2(i,j)=check1(i,:)*M(:,j);
        if(abs(check2(i,j)-saveA(i,j))>10^-10)
            disp("This inverse is incorrect")
        end
    end
end

    
        