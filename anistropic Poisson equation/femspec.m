%求解有限元方法线高斯迭代矩阵的谱半径
function r=femspec(n,eps,omega)
%omega为区域旋转角度
%线高斯的迭代矩阵依然是inv（A）*B，A为分块下三角矩阵，B为分块上三角阵并且取负号。以下分别计算A，B矩阵
n=2^n-1;%网格规模
c = cos(omega);s = sin(omega);
gamma = 2 + 2*c*s + eps*(2 - 2*c*s);
alpha = -(c*s + c*c) - eps*(s*s - c*s);
beta = -(s*s + c*s) - eps*(c*c - c*s);
delta = (1 - eps)*c*s;
A=sparse(gamma*diag(ones(n,1),0)+ alpha*diag(ones(n-1,1),1)+alpha*diag(ones(n-1,1),-1));
A=kron(eye(n),A);
Bt=sparse(beta*diag(ones(n,1),0)+delta*diag(ones(n-1,1),-1));
B=Bt';
A(n+1:n^2,1:n^2-n)=A(n+1:n^2,1:n^2-n)+kron(eye(n-1),Bt);
B=-1*kron(diag(ones(n-1,1),1),B);
x=ones(n^2,1);
count=0;
r0=norm(A\(B*x))/norm(x);
%以下使用幂法对inv（A）*B的谱半径进行求解
while 1
    x=A\(B*x);
    x=x/norm(x,'inf');
    count=count+1;
    if count==50
        r=norm(A\(B*x))/norm(x);
        if abs(r-r0)<10^-5
            break;
        else
            r0=r;
            count=0;
        end
    end
end



