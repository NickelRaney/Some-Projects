function y=r(v);
%本函数的目的是随机采样，在给定最优前景区域后进行分块采样。
%环绕区域采样可能会在另外一个函数里给出
n=length(v);
y=v;
flag=zeros(n);
i=1;
while v(i)~=0
    flag(v(i))=1;
    i=i+1;
end
%以上是录入当前最有前景区域，flag的目的是帮助随机
res=find(flag==0);
perm=randperm(n-i+1);
for j=i:n
    y(j)=res(perm(j-i+1));
end
end


