function [U,S,V]=svd_ex(A)
Err=realmax; %Preventing for stopping the program.
sizeA=size(A);
cntmax=100*max(sizeA);cnt=0; %For precision
U=eye(sizeA(1));S=A';V=eye(sizeA(2));
while cnt<cntmax
    [Q,S]=QR_H(S');U=U*Q; %QR_H is the decomposition using Gram_Schmidt
    [Q,S]=QR_H(S');V=V*Q;
    e=triu(S,1);
    E=norm(e(:));
    F=norm(diag(S));
    if F==0, F=1; end
    Err=E/F;
    cnt=cnt+1;
end
%For fixing the sign of Sigma
sigma=diag(S);
for i=1:length(sigma)
    sigmad=sigma(i);
    sigma(i,i)=abs(sigmad);
    if sigmad<0
        u(:,i)=-u(:,i);
    end
end
return;

