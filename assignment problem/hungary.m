%Hungary-algorithm
function [t,cost,rm,cm]=hungary(c)
tic;
n=length(c);
cost=0;
c0=c;

%Ԥ����
mcol=min(c,[],1);
c=c-ones(n,1)*mcol;
mrow=min(c,[],2);
c=c-mrow*ones(1,n);
%�Ƿ���n������Ԫ

%������Ҫһ��ѭ�����������ֱ����������Ԫ����Ϊn
while 1
    c2=c;
    rm=zeros(1,n);
    cm=zeros(1,n);
    [zcol,zrow]=find(c'==0);
    zrow0=zrow;
    rn=zeros(1,n);%ͳ��ÿ�е���Ԫ�ĸ���
    for i=1:n
        rn(i)=length(find(zrow==i));
    end
    count=1;
    while max(rn)~=0
        frn=rn;
        frn(rn==0)=max(rn(:))+1;
        [minx,ind]=min(frn);%��С����Ԫ�б�
        rm(count)=ind;
        v=find(zrow==ind);
        cm(count)=zcol(v(1));%���б��
        count=count+1;
        rn(ind)=0;
        v=find(zcol==zcol(v(1)));%�ҳ��ñ��Ԫ�����е���Ԫ
        for i=1:length(v)
            if zrow(v(i))~=ind
                rn(zrow(v(i)))=rn(zrow(v(i)))-1;
                zrow(v(i))=0;
            end
        end
    end
    
    flag=zeros(1,n);
    for i=1:n
        if rm(i)~=0
            flag(rm(i))=1;
        end
    end
    rindex=find(flag==0);
    
    %�����㷨�ó��˱���ǵ�Ԫ��
    if isempty(rindex)==1
        cost=0;
        result=sortrows([rm;cm]')';
        result
        for k=1:n
            c0(k,result(2,k))
            cost=cost+c0(k,result(2,k));
        end
        cost
        toc;
        t=toc;
        break;
    end
    
  
    zrow=zrow0;
    
    
    
    %���½��д�
    cindex=zeros(1,n);
    for i=1:length(rindex)
        index2=find(zrow==rindex(i));%����������Ԫ
        
        if isempty(index2)==0%���Ԫ�ش���
            for j=1:length(index2)
                if isempty(find(cindex==zcol(index2(j))))==1%���cindex���޸���
                    no=find(cindex==0);
                    cindex(no(1))=zcol(index2(j));
                end
            end
            index3=zeros(1,n);%��¼���λ�ȥ������
            index4=zeros(1,n);%��¼���λ�ȥ������
            vcol=zcol(index2);%��Ӧ����Ԫ����
            while 1
                count=0;%count�ж��Ƿ�򹴣����ж�ͣ�� 
                for j=1:length(vcol)
                    if isempty(find(cm==vcol(j)))==0%����������б��Ԫ�أ������һ��
                        r=rm(find(cm==vcol(j)));%�ñ��Ԫ��Ӧ������
                        if isempty(find(rindex==r))==1%��δ�����
                            no=length(rindex);
                            rindex0=zeros(1,no+1);
                            rindex0(1:no)=rindex;
                            rindex0(no+1)=r;
                            rindex=rindex0;
                            
                            count=count+1;
                            no=find(index3==0);
                            index3(no(1))=r;
                        end
                    end
                end
                for j=1:length(index3)
                    no=find(zrow==index3(j));%�����еķ���Ԫ������
                    for k=1:length(no)
                        if zcol(no(k))~=cm(find(rm==index3(j)))%��Ϊ���Ԫ
                            if isempty(find(cindex==zcol(no(k))))==1
                                noo=find(cindex==0);
                                cindex(noo(1))=zcol(no(k));
                                noo=find(index4==0);
                                index4(noo(1))=zcol(no(k));
                                count=count+1;
                            end
                        end
                    end
                end
                if count==0
                    break;
                end
                no=find(index4==0);
                vcol=index4(1,1:no(1)-1);
            end
        end
    end
   
    
    maxc=max(max(c));
    for i=1:n
        if isempty(find(rindex==i))==1
            c2(i,:)=maxc+1;
        end
    end
    no=find(cindex==0);
    if no(1)~=0  
        c2(:,cindex(1:no(1)-1))=maxc+1;
    end
    minc=min(min(c2));
    
    for i=1:n
        if isempty(find(rindex==i))==0
            c(i,:)=c(i,:)-minc;
        end
    end
    if no(1)~=0  
        c(:,cindex(1:no(1)-1))=c(:,cindex(1:no(1)-1))+minc;
    end
end


















