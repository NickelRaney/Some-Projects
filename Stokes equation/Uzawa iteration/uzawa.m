function [u,v,p]=uzawa(A,N,p,f,g,d,iter,w)
n=2^N;
h=1/n;
freal=f;
greal=g;
dreal=d;
%f-Bp
for count=1:iter
    f=freal-(p(2:n,:)-p(1:n-1,:))/h;
    g=greal-(p(:,2:n)-p(:,1:n-1))'/h;
    fr=f(:);
    gr=g(:);
    ur=A\fr;
    vr=A\gr;
    u=reshape(ur,n-1,n);
    v=reshape(vr,n-1,n);
    d=dreal;
    d(1:n-1,:)=d(1:n-1,:)+u(:,:)/h;
    d(2:n,:)=d(2:n,:)-u(:,:)/h;
    d(:,1:n-1)=d(:,1:n-1)+v(:,:)'/h;
    d(:,2:n)=d(:,2:n)-v(:,:)'/h;
    p=p-w*d;
end



