function y=r(v);
%��������Ŀ��������������ڸ�������ǰ���������зֿ������
%��������������ܻ�������һ�����������
n=length(v);
y=v;
flag=zeros(n);
i=1;
while v(i)~=0
    flag(v(i))=1;
    i=i+1;
end
%������¼�뵱ǰ����ǰ������flag��Ŀ���ǰ������
res=find(flag==0);
perm=randperm(n-i+1);
for j=i:n
    y(j)=res(perm(j-i+1));
end
end


