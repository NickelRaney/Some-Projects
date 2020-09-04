function f=ana_sol(X,Y,eps)
syms x y
s3=sqrt(3);
u=((s3*x+y-s3/2).*(s3*x+y-s3/2-2).*(x-s3*y-1/2).*(x-s3*y+3/2));
dfx=diff(u,x);
dfy=diff(u,y);
fun=-diff(dfx,x)-eps*diff(dfy,y);
myfun=matlabFunction(fun);
f=myfun(X,Y);
% mesh(X,Y,f)
% f=subs(fun,{x,y},{1,2})
end