function [H]=house(a)
[n,m]=size(a);
e1=[1;zeros(1,n-1)'];
k=-sign(a(1))*norm(a);
lambda=sqrt(2*norm(a)*(norm(a)+abs(a(1))));
v=(a-k*e1)./lambda;
H=eye(n)-2*v*v';
end