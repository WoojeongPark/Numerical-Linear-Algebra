N=input("What are the dimension of the matrix n and max_iter?")
max_iter=input("What are the dimension of the matrix n and max_iter?")
h=1/N; iter=0; tol=1.0e-15; A=zeros(N);
for j=1:N
    A(j,j)=2/h^2;
end
for j=1:N-1
    A(j,j+1)=-1/h^2;
    A(j+1,j)=-1/h^2;
end
A(N,1)=-1/h^2;A(1,N)=-1/h^2;
q=zeros(N,1); q(1)=1; %Initial approximated eigenvector q
prev_q=q;
%Using Inverse Iterative Method to find minimum eigenvalue and eigenvector
q=zeros(N,1); q(1)=1; %Initial approximated eigenvector q
prev_q=q;lambda=zeros(41,1);
for j=0:10:400
    T=A-(j+0.1)*eye(N);
    for iter=1:max_iter
        z=T\prev_q;
        q=z/norm(z);
        diff_v=norm(q-prev_q);
        prev_q=q;
        if(diff_v<tol)
            disp("Converged at iter=")
            iter
            break
        end
    end
    lambda(j/10+1)=q'*A*q;
end
disp("Computed eigenvalue is")
lambda
       