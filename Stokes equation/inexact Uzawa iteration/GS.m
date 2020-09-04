function u=GS(N,u,f)
%N���������u,v�Ĺ�ģΪ(2^n-1)*2^n,p�Ĺ�ģΪ2^n*2^n
n=2^N;
h=1/n;
%Gauss-Seidel
%����F-Uu
f(1:n-2,:)=f(1:n-2,:)+u(2:n-1,:)/h^2;
f(:,1:n-1)=f(:,1:n-1)+u(:,2:n)/h^2;
%����inv(D-L)*(F-Uu)
i=1;
j=1;
u(i,j)=h^2*f(i,j)/3;
for i=2:n-1
    u(i,j)=h^2*f(i,j)/3+u(i-1,j)/3;
end
for j=2:n-1
    i=1;
    u(i,j)=h^2*f(i,j)/4+u(i,j-1)/4;
    for i=2:n-1
        u(i,j)=h^2*f(i,j)/4+u(i,j-1)/4+u(i-1,j)/4;
    end
end
j=n;
i=1;
u(i,j)=h^2*f(i,j)/3+u(i,j-1)/3;
for i=2:n-1
    u(i,j)=h^2*f(i,j)/3+u(i,j-1)/3+u(i-1,j)/3;
end













