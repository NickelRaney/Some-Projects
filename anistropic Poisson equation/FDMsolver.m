function [u,k,resVec,iniTime,solTime] = FDMsolver(n,eps,w,F)
tic
u = zeros(n-1);

c = cos(w);
s = sin(w);
e = ones(n-1, 1);

T = spdiags([-e 2*e -e], -1: 1, n-1, n-1);
R = spdiags([-e e], [-1 1], n-1, n-1)/2;
Ah = (c^2+eps*s^2)*kron(speye(n-1), T) + (s^2+eps*c^2)*kron(T, speye(n-1)) + (1-eps)*2*c*s*kron(R, R);

% F = @(x,y)(-2*(c^4*y.^2+s^3*y.*(-1+s*y)+c*s*(1+3*s*x).*(-1+2*s*y)-... 
%    c^3*(y+6*s*x.*y)+c^2*s*(3*x+6*s*x.^2+2*y.*(1-2*s*y))))+eps*(...
%    -2*(c^4*x.^2+s^3*x.*(1+s*x)-c*s*(1+2*s*x).*(-1+3*s*y)+... 
%    c^3*x.*(-1+6*s*y)-c^2*s*(2*x+4*s*x.^2+3*y.*(1-2*s*y))));
% U = @(X1, Y1) (Y1*cos(w)-X1*sin(w)).*(-Y1*cos(w)+X1*sin(w)+1).*(X1*cos(w)+Y1*sin(w)).*(-X1*cos(w)-Y1*sin(w)+1);


A = (c^2+eps*s^2)*T + (s^2+eps*c^2)*2*speye(n-1);
B = -(s^2+eps*c^2)*speye(n-1) + (1-eps)*2*s*c*R/2;

x = 1/n: 1/n: 1-1/n;
y = 1/n: 1/n: 1-1/n;

[X, Y]=ndgrid(x, x);
f = F(c*X - s*Y, s*X + c*Y);
iniTime = toc;

r0 = norm(f, 'fro');
r1 = r0;
iter = 0;
resVec = zeros(80000,1);

tic
while r1/r0 > 1e-6
    resVec(iter+1) = r1;
    iter = iter+1;
    
    u(: , 1) = A\(f(: , 1)/n^2 - B*u(: , 2));
    
    for i = 2: n-2
        u(: , i) = A\(f(:, i)/n^2 - B'*u(: , i-1) - B*u(: , i+1));
    end
    
    u(: , n-1) = A\(f(:, n-1)/n^2 - B'*u(: , n-2));
    
    %     r1 = norm(f - ((c^2+eps*s^2)*T*u + (s^2+eps*c^2)*u*T + (1-eps)*2*c*s*R*u*R)*n^2, 'fro');
    r1 = norm(f(: ) - Ah*u(: )*n^2);
end
k = iter;
resVec(k+1) = r1;
resVec = resVec(1:k+1);
solTime = toc;

v = zeros(n+1);
v(2: n, 2: n) = u;
u = v;
