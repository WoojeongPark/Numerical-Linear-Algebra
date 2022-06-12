n=input("What are the dimension of the matrix n and max_iter?")
max_iter=input("What are the dimension of the matrix n and max_iter?")
h=1/n; iter=0; tol=1.0e-15; A=zeros(n-1);
for j=1:n-1
    A(j,j)=2/h^2;
end
for j=1:n-2
    A(j,j+1)=-1/h^2;
    A(j+1,j)=-1/h^2;
end
q=zeros(n-1,1); q(1)=1; %Initial approximated eigenvector q
prev_q=q;
%Using Power Method to find maximum eigenvalue and eigenvector
for iter=1:max_iter
    q=A*q;
    q=q/norm(q);
    lambda=q'*A*q;
    diff_v=norm(q-prev_q);
    prev_q=q;
    if(diff_v<tol)
        disp("Converged at iter=")
        iter
        break
    end
end
disp("Computed maximum eigenvalue is")
max_lambda=lambda
disp("Computed maximum eigenvector is")
max_v=q
if(A*max_v-max_lambda*max_v>tol) 
    disp("Something is wrong")
end
disp("-----------------------------------------------------------------")
%Using Inverse Iterative Method to find minimum eigenvalue and eigenvector
q=zeros(n-1,1); q(1)=1; %Initial approximated eigenvector q
prev_q=q;
for iter=1:max_iter
    z=A\prev_q;
    q=z/norm(z);
    diff_v=norm(q-prev_q);
    prev_q=q;
    if(diff_v<tol)
        disp("Converged at iter=")
        iter
        break
    end
end
lambda=q'*A*q;
disp("Computed minimum eigenvalue is")
min_lambda=lambda
disp("Computed minimum eigenvector is")
min_v=q
if(A*min_v-min_lambda*min_v>tol) 
    disp("Something is wrong")
end

       