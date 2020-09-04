%Nested-dissection algorithm
function [t,v,cost]=nested_dissection(c)
%c is cost matrix, t is running time, v is assignment, cost is total cost.
tic;
n=length(c);
v=zeros(1,n);
flag=zeros(1,n);
level=0;
while level~=n-1
    choice=find(flag==0);
    %sampling in the dissected optimal region
    val=zeros(1,n-level);
    mino=0;
    min=n;
    for i=1:n-level 
        tot=0;
        v(level+1)=choice(i);
        mintot=n;
        for j=1:ceil((n-level)^1.5)%number of samples
           sample=r(v);
           for k=1:n
               tot=tot+c(k,sample(k));
           end
           if tot<mintot
               mintot=tot;
           end
        end
        val(i)=mintot;%find the minimum of samples
        if min>val(i)
            mino=i;
            min=val(i);
        end
    end
    %sampling in surrounding region
    if level==0
        v(level+1)=choice(mino);
        level=level+1;
        flag(choice(mino))=1;
    else
        tot=0;
        mintot=n;
        for i=1:ceil((n-level)^1.5)%number of sampling
            sample=r(zeros(1,n));
            while sample(1:level)==v(1:level)
                sample=r(zeros(1,n));
            end
            for k=1:n
               tot=tot+c(k,sample(k));
            end
            if tot<mintot
                mintot=tot;
            end
        end
        valo=mintot;
        if valo<val(mino)
            flag(v(level))=0;
            v(level)=0;
            level=level-1;
        else
            v(level+1)=choice(mino);
            level=level+1;
            flag(choice(mino))=1;
        end
    end
end
v=r(v);
cost=0;
for k=1:n
    cost=cost+c(k,v(k));
end
toc;
t=toc;


        