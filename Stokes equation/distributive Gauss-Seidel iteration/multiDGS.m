function [err,t,count]=multiDGS(iter1,iter2,l,n)
%DGS多重网格,iter1为下降磨光次数，iter2为上升磨光次数，l为下降层数，最大为n-2，网格规模为2^n
u=cell(n,1);
v=cell(n,1);
f=cell(n,1);
g=cell(n,1);
p=cell(n,1);
d=cell(n,1);

for i=n:-1:n-l
    u{i}=zeros(2^i-1,2^i);
    v{i}=zeros(2^i-1,2^i);
    f{i}=zeros(2^i-1,2^i);
    g{i}=zeros(2^i-1,2^i);
    p{i}=zeros(2^i,2^i);
    d{i}=zeros(2^i,2^i);
end
h=1/2^n;
%f，g为方程右端项
for j=1:2^n
    for i=1:2^n-1
        f{n}(i,j)=funcf(i*h,(j-0.5)*h);
        g{n}(i,j)=funcg((j-0.5)*h,i*h);
    end
end
%f，g加上边界条件
for i=1:2^n-1
    f{n}(i,1)=f{n}(i,1)+b(i*h)/h;
    f{n}(i,2^n)=f{n}(i,2^n)-b(i*h)/h;
    g{n}(i,1)=g{n}(i,1)-b(i*h)/h;
    g{n}(i,2^n)=g{n}(i,2^n)+b(i*h)/h;
end
fr=f{n};
gr=g{n};
%求初始误差r0
r0=sqrt(norm(fr,'fro')^2+norm(gr,'fro')^2);
%下面开始多重网格
%下降
tic;
count=0;
while 1
    for i=n:-1:n-l+1
        if i~=n
            u{i}=zeros(2^i-1,2^i);
            v{i}=zeros(2^i-1,2^i);
            p{i}=zeros(2^i,2^i);
        end
        for j=1:iter1
            [u{i},v{i},p{i}]=DGS(i,u{i},v{i},p{i},f{i},g{i},d{i});
        end
        [fr,gr,dr]=res(i,u{i},v{i},p{i},f{i},g{i},d{i});
        %限制
        f{i-1}=0.125*fr(1:2:2^i-3,1:2:2^i-1)+0.125*fr(1:2:2^i-3,2:2:2^i)+0.125*fr(3:2:2^i-1,1:2:2^i-1)+0.125*fr(3:2:2^i-1,2:2:2^i)+0.25*fr(2:2:2^i-2,1:2:2^i-1)+0.25*fr(2:2:2^i-2,2:2:2^i);
        g{i-1}=0.125*gr(1:2:2^i-3,1:2:2^i-1)+0.125*gr(1:2:2^i-3,2:2:2^i)+0.125*gr(3:2:2^i-1,1:2:2^i-1)+0.125*gr(3:2:2^i-1,2:2:2^i)+0.25*gr(2:2:2^i-2,1:2:2^i-1)+0.25*gr(2:2:2^i-2,2:2:2^i);
        d{i-1}=0.25*dr(1:2:2^i-1,1:2:2^i-1)+0.25*dr(1:2:2^i-1,2:2:2^i)+0.25*dr(2:2:2^i,1:2:2^i-1)+0.25*dr(2:2:2^i,2:2:2^i);
    end
    i=n-l;
    u{i}=zeros(2^i-1,2^i);
    v{i}=zeros(2^i-1,2^i);
    p{i}=zeros(2^i,2^i);
    for j=1:iter1
        [u{i},v{i},p{i}]=DGS(i,u{i},v{i},p{i},f{i},g{i},d{i});
    end
    %上升
    for i=n-l+1:n
        u{i}(1:2:2^i-3,1:2:2^i-1)=u{i}(1:2:2^i-3,1:2:2^i-1)+0.5*u{i-1};
        u{i}(1:2:2^i-3,2:2:2^i)=u{i}(1:2:2^i-3,2:2:2^i)+0.5*u{i-1};
        u{i}(1:2:2^i-3,2:2:2^i)=u{i}(1:2:2^i-3,2:2:2^i)+0.5*u{i-1};
        u{i}(3:2:2^i-1,2:2:2^i)=u{i}(3:2:2^i-1,2:2:2^i)+0.5*u{i-1};
        u{i}(2:2:2^i-2,1:2:2^i-1)=u{i}(2:2:2^i-2,1:2:2^i-1)+u{i-1};
        u{i}(2:2:2^i-2,2:2:2^i)=u{i}(2:2:2^i-2,2:2:2^i)+u{i-1};
        
        v{i}(1:2:2^i-3,1:2:2^i-1)=v{i}(1:2:2^i-3,1:2:2^i-1)+0.5*v{i-1};
        v{i}(1:2:2^i-3,2:2:2^i)=v{i}(1:2:2^i-3,2:2:2^i)+0.5*v{i-1};
        v{i}(1:2:2^i-3,2:2:2^i)=v{i}(1:2:2^i-3,2:2:2^i)+0.5*v{i-1};
        v{i}(3:2:2^i-1,2:2:2^i)=v{i}(3:2:2^i-1,2:2:2^i)+0.5*v{i-1};
        v{i}(2:2:2^i-2,1:2:2^i-1)=v{i}(2:2:2^i-2,1:2:2^i-1)+v{i-1};
        v{i}(2:2:2^i-2,2:2:2^i)=v{i}(2:2:2^i-2,2:2:2^i)+v{i-1};
        
        p{i}(1:2:2^i-1,1:2:2^i-1)=p{i}(1:2:2^i-1,1:2:2^i-1)+p{i-1};
        p{i}(1:2:2^i-1,2:2:2^i)=p{i}(1:2:2^i-1,2:2:2^i)+p{i-1};
        p{i}(2:2:2^i,1:2:2^i-1)=p{i}(2:2:2^i,1:2:2^i-1)+p{i-1};
        p{i}(2:2:2^i,2:2:2^i)=p{i}(2:2:2^i,2:2:2^i)+p{i-1};
        for j=1:iter2
            [u{i},v{i},p{i}]=DGS(i,u{i},v{i},p{i},f{i},g{i},d{i});
        end
    end
    %res函数求F-Bp-Au
    [fr,gr,dr]=res(n,u{n},v{n},p{n},f{n},g{n},d{n});
    %r3为一个cycle之后的误差
    r3=sqrt(norm(fr,'fro')^2+norm(gr,'fro')^2+norm(dr,'fro')^2);
    if r3<10^-8*r0
        break;
    end
    count=count+1;
end
%urr,vrr,prr为真解
t=toc;
urr=zeros(2^n-1,2^n);
vrr=zeros(2^n-1,2^n);
prr=zeros(2^n,2^n);
h=1/2^n;
for j=1:2^n
    for i=1:2^n-1
        urr(i,j)=ur(i*h,(j-0.5)*h);
        vrr(i,j)=vr((j-0.5)*h,i*h);
    end
end
%求数值解与真解误差
err=h*sqrt(norm(urr-u{n},'fro')^2+norm(vrr-v{n},'fro')^2);









