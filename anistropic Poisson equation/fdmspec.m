%求解九点差分格式线高斯迭代矩阵谱半径
function r=fdmspec(n,eps,omega)
%omega为区域旋转角度
%线高斯的迭代矩阵依然是inv（A）*B，A为分块下三角矩阵，B为分块上三角阵并且取负号。以下分别计算A，B矩阵
n=2^n-1;%网格规模
c = cos(omega);s = sin(omega);
e = ones(n, 1);
T = spdiags([-e 2*e -e], -1: 1, n, n);
R = spdiags([-e e], [-1 1], n, n)/2;
A = (c^2+eps*s^2)*T + (s^2+eps*c^2)*2*speye(n);
B = -(s^2+eps*c^2)*speye(n) + (1-eps)*2*s*c*R/2;
A=kron(eye(n),A);
A(n+1:n^2,1:n^2-n)=A(n+1:n^2,1:n^2-n)+kron(eye(n-1),B');
B=-1*kron(diag(ones(n-1,1),1),B);
x=ones(n^2,1);
count=0;
%以下使用幂法对inv（A）*B的谱半径进行求解
r0=norm(A\(B*x))/norm(x);
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




