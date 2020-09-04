function [u,k,resVec,iniTime,solTime] = FDM45solver(N,eps,F)
tic;
w = pi/4;
c = cos(w); s = sin(w);   

f = zeros(2*N+1);     % 右端项
u = zeros(2*N+1);     % 数值解

%% 生成每层矩阵
Ai = cell(N+1,1); 
for i = 2:N+1
    M = 2*i-3; 
    Ai{i} = sparse(2*(1+eps)*diag(ones(M,1))-diag(ones(M-1,1),1)-diag(ones(M-1,1),-1)); 
end

%% 生成右端项
for i = 2:N+1
    y = 2+1/N-i/N;
    for j = N+3-i:N-1+i
        x = -1-1/N+j/N;
        f(i,j) = 1/2/N^2*F(x/sqrt(2),y/sqrt(2));
%         f(i,j) = 1/N^2;
    end
end

for i = N:-1:2
    y = -1/N+i/N;
    for j = N+3-i:N-1+i
        x = -1-1/N+j/N;
        f(2*N+2-i,j) = 1/2/N^2*F(x/sqrt(2),y/sqrt(2));
%         f(2*N+2-i,j) = 1/N^2;
    end
end

iniTime = toc;

%% 生成初始残差
res0 = 0; 
for i = 2:N+1
    res0 = res0+norm(f(i,N+3-i:N-1+i)'-Ai{i}*u(i,N+3-i:N-1+i)'...
        +eps*(u(i-1,N+3-i:N-1+i)'+u(i+1,N+3-i:N-1+i)'),2)^2;
end
for i = N:-1:2
    res0 = res0+norm(f(2*N+2-i,N+3-i:N-1+i)'-Ai{i}*u(2*N+2-i,N+3-i:N-1+i)'...
        +eps*(u(2*N+1-i,N+3-i:N-1+i)'+u(2*N+3-i,N+3-i:N-1+i)'),2)^2;
end
res0 = sqrt(res0);
res = res0;

%% FDM45求解器
tic;
k = 0;
resVec = zeros(100000,1);
while res > 1e-6*res0
    resVec(k+1) = res;
    
    res = 0;
    for i = 2:N+1
        vec = Ai{i}\(f(i,N+3-i:N-1+i)'+eps*(u(i-1,N+3-i:N-1+i)'+u(i+1,N+3-i:N-1+i)')); 
        u(i,N+3-i:N-1+i) = vec'; 
    end
    for i = N:-1:2
        vec = Ai{i}\(f(2*N+2-i,N+3-i:N-1+i)'+eps*(u(2*N+1-i,N+3-i:N-1+i)'+u(2*N+3-i,N+3-i:N-1+i)')); 
        u(2*N+2-i,N+3-i:N-1+i) = vec';
    end
    
    for i = 2:N+1
        res = res+norm(f(i,N+3-i:N-1+i)'-Ai{i}*u(i,N+3-i:N-1+i)'...
            +eps*(u(i-1,N+3-i:N-1+i)'+u(i+1,N+3-i:N-1+i)'),2)^2;
    end
    for i = N:-1:2
        res = res+norm(f(2*N+2-i,N+3-i:N-1+i)'-Ai{i}*u(2*N+2-i,N+3-i:N-1+i)'...
            +eps*(u(2*N+1-i,N+3-i:N-1+i)'+u(2*N+3-i,N+3-i:N-1+i)'),2)^2;
    end
    
    res = sqrt(res);
    k = k+1;
end
resVec(k+1) = res;
resVec = resVec(1:k+1);
solTime = toc;


end