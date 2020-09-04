%�������Ԫ�����߸�˹����������װ뾶
function r=femspec(n,eps,omega)
%omegaΪ������ת�Ƕ�
%�߸�˹�ĵ���������Ȼ��inv��A��*B��AΪ�ֿ������Ǿ���BΪ�ֿ�����������ȡ���š����·ֱ����A��B����
n=2^n-1;%�����ģ
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
%����ʹ���ݷ���inv��A��*B���װ뾶�������
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



