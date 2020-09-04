function [u,v,p]=ineuzawa(N,p,f,g,d,iter,w,iter1)
n=2^N;
h=1/n;
freal=f;
greal=g;
dreal=d;
%f-Bp
for count=1:iter
    f=freal-(p(2:n,:)-p(1:n-1,:))/h;
    g=greal-(p(:,2:n)-p(:,1:n-1))'/h;
    if N>2
    u=mg(f,iter1,iter1,N-2,N);
    v=mg(g,iter1,iter1,N-2,N);
    else
        u=zeros(3,4);
        v=zeros(3,4);
        u=GS(2,u,f);
        v=GS(2,u,f);
    end
    d=dreal;
    d(1:n-1,:)=d(1:n-1,:)+u(:,:)/h;
    d(2:n,:)=d(2:n,:)-u(:,:)/h;
    d(:,1:n-1)=d(:,1:n-1)+v(:,:)'/h;
    d(:,2:n)=d(:,2:n)-v(:,:)'/h;
    p=p-w*d;
end



