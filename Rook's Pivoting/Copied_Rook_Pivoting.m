A=input('LU분해를 알고 싶은 정사각행렬을 입력하시오: ');
[n,m]=size(A);
if(n~=m)
    disp('정사각행렬이 아닙니다.');
    return;
end
L=zeros(n,n); P_row=1:n; P_col=1:n;
r=1; c=1; r0=0; c0=0;
for i=1:n
    while r~=r0 || c~=c0
        r0=r;
        c0=c;
        [rmax,c]=max(abs(A(r+i-1,i:n)));
        [cmax,r]=max(abs(A(i:n,c+i-1)));
    end
    %가장 큰 pivot을 찾았으니, 배치한다.
    r=r+i-1;
    c=c+i-1;
    temp=P_row(i);
    P_row(i)=P_row(r);
    P_row(r)=temp;
    temp=P_col(i);
    P_col(i)=P_col(c);
    P_col(c)=temp;
    for j=1:n
        temp=A(i,j);
        A(i,j)=A(r,j);
        A(r,j)=temp;
    end
    for k=1:n
        temp=A(k,i);
        A(k,i)=A(k,c);
        A(k,c)=temp;
    end
    if A(i,i)~=0
        for k=i+1:n
            A(k,i)=A(k,i)/A(i,i);
        end
        for j=i+1:n
            for k=i+1:n
                A(k,j)=A(k,j)-A(k,i)*A(i,j);
            end
        end
    end
end
disp("결과물")
P
A