function [err,t,count]=multiuzawa(iter1,iter2,l,n,w)
%Uzawa多重网格,iter1为下降磨光次数，iter2为上升磨光次数，l为下降层数，最大为n-2，网格规模为2^n
u=cell(n,1);
v=cell(n,1);
f=cell(n,1);
g=cell(n,1);
p=cell(n,1);
d=cell(n,1);
A=cell(n,1);
for le=2:n
    N=2^le;
    h=1/N;
    ex=ones(N-1,1);
    ey=ones(N,1);
    tx=spdiags([-ex 2*ex -ex],-1:1,N-1,N-1);
    ty=spdiags([-ey [0;ones(N-2,1);0]+ey -ey],-1:1,N,N);
    A{le}=kron(speye(N),tx)/h^2+kron(ty,speye(N-1))/h^2;
end
for i=n:-1:n-l
    u{i}=zeros(2^i-1,2^i);
    v{i}=zeros(2^i-1,2^i);
    f{i}=zeros(2^i-1,2^i);
    g{i}=zeros(2^i-1,2^i);
    p{i}=zeros(2^i,2^i);
    d{i}=zeros(2^i,2^i);
end
h=1/2^n;
for j=1:2^n
    for i=1:2^n-1
        f{n}(i,j)=funcf(i*h,(j-0.5)*h);
        g{n}(i,j)=funcg((j-0.5)*h,i*h);
    end
end
for i=1:2^n-1
    f{n}(i,1)=f{n}(i,1)+b(i*h)/h;
    f{n}(i,2^n)=f{n}(i,2^n)-b(i*h)/h;
    g{n}(i,1)=g{n}(i,1)-b(i*h)/h;
    g{n}(i,2^n)=g{n}(i,2^n)+b(i*h)/h;
end
fr=f{n};
gr=g{n};
r0=sqrt(norm(fr,'fro')^2+norm(gr,'fro')^2);
tic;
count=0;
while 1
    for i=n:-1:n-l+1
        if i~=n
            u{i}=zeros(2^i-1,2^i);
            v{i}=zeros(2^i-1,2^i);
            p{i}=zeros(2^i,2^i);
        end
        [u{i},v{i},p{i}]=uzawa(A{i},i,p{i},f{i},g{i},d{i},iter1,w);
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
    [u{i},v{i},p{i}]=uzawa(A{i},i,p{i},f{i},g{i},d{i},iter1,w);
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
        [u{i},v{i},p{i}]=uzawa(A{i},i,p{i},f{i},g{i},d{i},iter2,w);
    end
    [fr,gr,dr]=res(n,u{n},v{n},p{n},f{n},g{n},d{n});
    r3=sqrt(norm(fr,'fro')^2+norm(gr,'fro')^2+norm(dr,'fro')^2);
    if r3<10^-8*r0
        break;
    end
    count=count+1;
end
count=count+1;
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
u{n}=u{n}-urr;
v{n}=v{n}-vrr;
err=h*sqrt(norm(u{n},'fro')^2+norm(v{n},'fro')^2);





