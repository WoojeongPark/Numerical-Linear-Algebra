tol=1.0e-13;
for n=4:7
    A=rand(n);
    [u,s,v]=svd_ex(A);
    tol=tol/10;
    if(norm(A-u*s*v')>tol)
    fprintf('The computed solution seems to be wrong at %f \n',n);
    end
end
