function [u,v,p]=DGS(N,u,v,p,f,g,D)
%N代表层数，u,v的规模为(2^n-1)*2^n,p的规模为2^n*2^n
n=2^N;
h=1/n;
%Gauss-Seidel

f=f+(p(1:n-1,:)-p(2:n,:))/h;
g=g+(p(:,1:n-1)-p(:,2:n))'/h;
f(1:n-2,:)=f(1:n-2,:)+u(2:n-1,:)/h^2;
g(1:n-2,:)=g(1:n-2,:)+v(2:n-1,:)/h^2;
f(:,1:n-1)=f(:,1:n-1)+u(:,2:n)/h^2;
g(:,1:n-1)=g(:,1:n-1)+v(:,2:n)/h^2;
%计算inv(D-L)*(F-Bp+Uu)
i=1;
j=1;
u(i,j)=h^2*f(i,j)/3;
v(i,j)=h^2*g(i,j)/3;
for i=2:n-1
    u(i,j)=h^2*f(i,j)/3+u(i-1,j)/3;
    v(i,j)=h^2*g(i,j)/3+v(i-1,j)/3;
end
for j=2:n-1
    i=1;
    u(i,j)=h^2*f(i,j)/4+u(i,j-1)/4;
    v(i,j)=h^2*g(i,j)/4+v(i,j-1)/4;
    for i=2:n-1
        u(i,j)=h^2*f(i,j)/4+u(i,j-1)/4+u(i-1,j)/4;
        v(i,j)=h^2*g(i,j)/4+v(i,j-1)/4+v(i-1,j)/4;
    end
end
j=n;
i=1;
u(i,j)=h^2*f(i,j)/3+u(i,j-1)/3;
v(i,j)=h^2*g(i,j)/3+v(i,j-1)/3;
for i=2:n-1
    u(i,j)=h^2*f(i,j)/3+u(i,j-1)/3+u(i-1,j)/3;
    v(i,j)=h^2*g(i,j)/3+v(i,j-1)/3+v(i-1,j)/3;
end


%更新速度和压力
%内部节点


for i=2:n-1
    for j=2:n-1
        r=-(u(i,j)-u(i-1,j))-(v(j,i)-v(j-1,i))-h*D(i,j);
        u(i,j)=u(i,j)+r/4;
        u(i-1,j)=u(i-1,j)-r/4;
        v(j,i)=v(j,i)+r/4;
        v(j-1,i)=v(j-1,i)-r/4;
        r=r/h;
        p(i,j)=p(i,j)+r;
        p(i+1,j)=p(i+1,j)-r/4;
        p(i-1,j)=p(i-1,j)-r/4;
        p(i,j+1)=p(i,j+1)-r/4;
        p(i,j-1)=p(i,j-1)-r/4;
    end
end
%下左右上更新四条边
j=1;
for i=2:n-1
    r=-(u(i,j)-u(i-1,j))-v(j,i)-h*D(i,j);
    u(i,j)=u(i,j)+r/3;
    u(i-1,j)=u(i-1,j)-r/3;
    v(j,i)=v(j,i)+r/3;
    r=r/h;
    p(i,j)=p(i,j)+4/3*r;
    p(i+1,j)=p(i+1,j)-r/3;
    p(i-1,j)=p(i-1,j)-r/3;
    p(i,j+1)=p(i,j+1)-r/3;
end
i=1;
for j=2:n-1
    r=-u(i,j)-(v(j,i)-v(j-1,i))-h*D(i,j);
    u(i,j)=u(i,j)+r/3;
    v(j,i)=v(j,i)+r/3;
    v(j-1,i)=v(j-1,i)-r/3;
    r=r/h;
    p(i,j)=p(i,j)+4/3*r;
    p(i+1,j)=p(i+1,j)-r/3;
    p(i,j-1)=p(i,j-1)-r/3;
    p(i,j+1)=p(i,j+1)-r/3;
end
i=n;
for j=2:n-1
    r=u(i-1,j)-(v(j,i)-v(j-1,i))-h*D(i,j);
    u(i-1,j)=u(i-1,j)-r/3;
    v(j,i)=v(j,i)+r/3;
    v(j-1,i)=v(j-1,i)-r/3;
    r=r/h;
    p(i,j)=p(i,j)+4/3*r;
    p(i-1,j)=p(i-1,j)-r/3;
    p(i,j-1)=p(i,j-1)-r/3;
    p(i,j+1)=p(i,j+1)-r/3;
end
j=n;
for i=2:n-1
    r=-(u(i,j)-u(i-1,j))+v(j-1,i)-h*D(i,j);
    u(i,j)=u(i,j)+r/3;
    u(i-1,j)=u(i-1,j)-r/3;
    v(j-1,i)=v(j-1,i)-r/3;
    r=r/h;
    p(i,j)=p(i,j)+4/3*r;
    p(i+1,j)=p(i+1,j)-r/3;
    p(i-1,j)=p(i-1,j)-r/3;
    p(i,j-1)=p(i,j-1)-r/3;
end
%更新四个角
%左下
i=1;j=1;
r=-u(i,j)-v(j,i)-h*D(i,j);
u(i,j)=u(i,j)+r/2;
v(j,i)=v(j,i)+r/2;
r=r/h;
p(i,j)=p(i,j)+2*r;
p(i+1,j)=p(i+1,j)-r/2;
p(i,j+1)=p(i,j+1)-r/2;
%右下
i=n;j=1;
r=u(i-1,j)-v(j,i)-h*D(i,j);
u(i-1,j)=u(i-1,j)-r/2;
v(j,i)=v(j,i)+r/2;
r=r/h;
p(i,j)=p(i,j)+2*r;
p(i-1,j)=p(i-1,j)-r/2;
p(i,j+1)=p(i,j+1)-r/2;
%左上
i=1;j=n;
r=-u(i,j)+v(j-1,i)-h*D(i,j);
u(i,j)=u(i,j)+r/2;
v(j-1,i)=v(j-1,i)-r/2;
r=r/h;
p(i,j)=p(i,j)+2*r;
p(i+1,j)=p(i+1,j)-r/2;
p(i,j-1)=p(i,j-1)-r/2;
%右上
i=n;j=n;
r=u(i-1,j)+v(j-1,i)-h*D(i,j);
u(i-1,j)=u(i-1,j)+r/2;
v(j-1,i)=v(j-1,i)+r/2;
r=r/h;
p(i,j)=p(i,j)+2*r;
p(i-1,j)=p(i-1,j)-r/2;
p(i,j-1)=p(i,j-1)-r/2;

end









