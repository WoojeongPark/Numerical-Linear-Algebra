max_iter=input("what is the max_iter?");
A=zeros(3);
A(1,1)=-261;A(1,2)=209;A(1,3)=-49;A(2,1)=-530;A(2,2)=422;A(2,3)=-98;A(3,1)=-800;A(3,2)=631;A(3,3)=-144;
v=[1 0 0]';prev_v=v;
tol=1.0e-15;
disp("Initial approximation of eigenvector is")
disp(prev_v)
for iter=1:max_iter
    z=A*v;
    v=z/norm(z);
    lambda=v'*A*v;
    diff_v=norm(v-prev_v);
    prev_v=v;
    if(diff_v<tol)
        disp("Converged at ")
        iter
        break
    end
end
disp("The maximum eigenvalue and eigenvector is")
lambda
v
if(A*v-lambda*v>tol) 
    disp("Something is wrong")
end