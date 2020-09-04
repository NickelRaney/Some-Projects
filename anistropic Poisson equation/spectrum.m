%这个程序求解N=64时三个特殊角度，不同eps大小的迭代矩阵并作图
x=zeros(9,1);
r1=zeros(9,1);
r2=zeros(9,1);
r3=zeros(9,1);
for i=1:9
    x(i)=i-1;
    eps=10^-(i-1);
    r1(i)=-1/log(femspec(6,eps,0));%同样可以使用命令fdmspec
    r2(i)=-1/log(femspec(6,eps,pi/6));
    r3(i)=-1/log(femspec(6,eps,pi/4));
end
plot(x,r1);
hold on;
plot(x,r2);
hold on;
plot(x,r3);
grid on;
title('N=64时谱半径与\epsilon关系')
xlabel('-log10(\epsilon)')
ylabel('-1/ln(spectrum)')
legend('\omega=0','\omega=\pi/6','\omega=\pi/4')
