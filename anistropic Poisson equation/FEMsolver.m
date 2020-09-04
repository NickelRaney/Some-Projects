function [u,k,resVec,iniTime,solTime] = FEMsolver(N,eps,w,F)
tic;
c = cos(w); s = sin(w);   

f = zeros(N+1);     % 右端项
u = zeros(N+1);     % 数值解

%% 构造右端项
for i = 1:N+1
    for j = 1:N+1
        x = (i-1)/N;y = (j-1)/N;
        x1 = c*x-s*y; y1 = s*x+c*y;
        f(i,j) = F(x1,y1)/N^2;
    end
end

%% 构造刚度矩阵
M = assembleMatrix(N,eps,w);

%% 计算矩阵的真解
% v(2:N,2:N) = reshape(M\reshape(f(2:N,2:N),(N-1)^2,1),N-1,N-1);

%% 计算刚度矩阵的子矩阵
A = ((2+2*c*s)+eps*(2-2*c*s))*diag(sparse(ones(N-1,1)))+...
    (-c*s-c^2-eps*(s^2-c*s))*(diag(sparse(ones(N-2,1)),1)+diag(sparse(ones(N-2,1)),-1));
B = (-s^2-c*s-eps*(c^2-c*s))*diag(sparse(ones(N-1,1)))+...
    (1-eps)*c*s*(diag(sparse(ones(N-2,1)),1));

iniTime = toc;

%% FEM求解器
tic;
k = 0;
res = norm(reshape(f(2:N,2:N),(N-1)^2,1)-M*reshape(u(2:N,2:N),(N-1)^2,1),2);
tol = 1e-6*res;
resVec = zeros(80000,1);
while res > tol
    resVec(k+1) = res;
%     fprintf('%d: %e\n',k,res);    % 输出每一步残差
    for j = 2:N
        u(2:N,j) = A\(f(2:N,j)-B'*u(2:N,j-1)-B*u(2:N,j+1));
    end
    k = k+1;
    res = norm(reshape(f(2:N,2:N),(N-1)^2,1)-M*reshape(u(2:N,2:N),(N-1)^2,1),2);
end
resVec(k+1) = res;
resVec = resVec(1:k+1);
solTime = toc;



end