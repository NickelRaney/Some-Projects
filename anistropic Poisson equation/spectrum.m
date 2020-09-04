%����������N=64ʱ��������Ƕȣ���ͬeps��С�ĵ���������ͼ
x=zeros(9,1);
r1=zeros(9,1);
r2=zeros(9,1);
r3=zeros(9,1);
for i=1:9
    x(i)=i-1;
    eps=10^-(i-1);
    r1(i)=-1/log(femspec(6,eps,0));%ͬ������ʹ������fdmspec
    r2(i)=-1/log(femspec(6,eps,pi/6));
    r3(i)=-1/log(femspec(6,eps,pi/4));
end
plot(x,r1);
hold on;
plot(x,r2);
hold on;
plot(x,r3);
grid on;
title('N=64ʱ�װ뾶��\epsilon��ϵ')
xlabel('-log10(\epsilon)')
ylabel('-1/ln(spectrum)')
legend('\omega=0','\omega=\pi/6','\omega=\pi/4')
