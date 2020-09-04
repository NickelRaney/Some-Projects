%数值实验主程序
n=5:8;

N=2.^n;
e=[1,1e-2,1e-4,1e-6,1e-8];

n=length(e);
m=length(N);

tT=zeros(n,m);
cnt=zeros(n,m);
err=zeros(n,m);
res=zeros(n,m);


for j=n:-1:1
    for i=1:m
        
        [tT(j,i),cnt(j,i),err(j,i),res(j,i)]=main_FDM_30(N(i),e(j));
        str=sprintf('N=%d, eps=%g. cnt=%d err=%g. using %g s. res=%g',N(i),e(j),cnt(j,i),err(j,i),tT(j,i),res(j,i));
        disp(str)
    end
end
save('FDM_FINAL_YEAH.mat')
plotData;%画图