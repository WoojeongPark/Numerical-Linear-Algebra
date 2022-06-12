max_iter=input("what is the max_iter?");
A=zeros(3);
A(1,1)=-261;A(1,2)=209;A(1,3)=-49;A(2,1)=-530;A(2,2)=422;A(2,3)=-98;A(3,1)=-800;A(3,2)=631;A(3,3)=-144;
v=[1 0 0]';prev_v=v;prev_A=A;
tol=1.0e-15;
disp("Initial approximation of eigenvector is")
disp(prev_v)
%Find the maximum eigenvalue and eigenvector : Exercise 4.11
%Find the eigenvector corresponding to the smallest eigenvalue
T=prev_A;
for iter=1:max_iter
    z=T\prev_v;
    v=z/norm(z);
    diff_v=norm(v-prev_v);
    prev_v=v;
    if(diff_v<tol)
        disp("Converged at ")
        iter
        break
    end
end
lambda=v'*A*v;
disp("The minimum eigenvalue and eigenvector is")
min_lambda=lambda
min_v=v
if(A*min_v-min_lambda*min_v>tol) 
    disp("Something is wrong")
end
%Find other eigenvector and eigenvalue.
v=[1 0 0]';prev_v=v;iter=0;
T=prev_A-6.5*eye(3);
%We already know max and min of eigenvalues, which are 10 and 3, so I will use the mean(10,3)=6.5
for iter=1:max_iter
    z=T\prev_v;
    v=z/norm(z);
    diff_v=norm(v-prev_v);
    prev_v=v;
    if(diff_v<tol)
        disp("Converged at ")
        iter
        break
    end
end
lambda=v'*A*v;
disp("The other eigenvalue and eigenvector is")
lambda
v
if(A*v-lambda*v>tol) 
    disp("Something is wrong")
end