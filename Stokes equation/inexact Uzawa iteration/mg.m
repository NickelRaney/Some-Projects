
function uu=mg(fr,iter1,iter2,l,n)
h=1/2^n;
u=cell(n,1);
f=cell(n,1);
for i=n:-1:n-l
    f{i}=zeros(2^i-1,2^i);
    u{i}=zeros(2^i-1,2^i);
end
f{n}=fr;
r0=norm(fr,'fro');
%下面开始多重网格
%下降
count=0;
while 1
    for i=n:-1:n-l+1
        if i~=n
            u{i}=zeros(2^i-1,2^i);
        end
        for j=1:iter1
            u{i}=GS(i,u{i},f{i});
        end
        fr=ires(i,u{i},f{i});
        %限制
        f{i-1}=0.125*fr(1:2:2^i-3,1:2:2^i-1)+0.125*fr(1:2:2^i-3,2:2:2^i)+0.125*fr(3:2:2^i-1,1:2:2^i-1)+0.125*fr(3:2:2^i-1,2:2:2^i)+0.25*fr(2:2:2^i-2,1:2:2^i-1)+0.25*fr(2:2:2^i-2,2:2:2^i);
    end
    i=n-l;
    u{i}=zeros(2^i-1,2^i);
    for j=1:iter1
        u{i}=GS(i,u{i},f{i});
    end
    %上升
    for i=n-l+1:n
        %细网格
        u{i}(1:2:2^i-3,1:2:2^i-1)=u{i}(1:2:2^i-3,1:2:2^i-1)+0.5*u{i-1};
        u{i}(1:2:2^i-3,2:2:2^i)=u{i}(1:2:2^i-3,2:2:2^i)+0.5*u{i-1};
        u{i}(1:2:2^i-3,2:2:2^i)=u{i}(1:2:2^i-3,2:2:2^i)+0.5*u{i-1};
        u{i}(3:2:2^i-1,2:2:2^i)=u{i}(3:2:2^i-1,2:2:2^i)+0.5*u{i-1};
        u{i}(2:2:2^i-2,1:2:2^i-1)=u{i}(2:2:2^i-2,1:2:2^i-1)+u{i-1};
        u{i}(2:2:2^i-2,2:2:2^i)=u{i}(2:2:2^i-2,2:2:2^i)+u{i-1};
        for j=1:iter2
            u{i}=GS(i,u{i},f{i});
        end
    end
    fr=ires(n,u{n},f{n});
    r3=norm(fr,'fro');
    urr=zeros(2^n+1,2^n);
    urr(2:2^n,:)=u{n};
    drr=(urr(1:2^n,:)-urr(2:2^n+1,:))*2^n;
    r4=norm(drr,'fro');
    if r3<0.000001*r4%在此处调节\tau
        break;
    end
    count=count+1;
end
uu=u{n};





















