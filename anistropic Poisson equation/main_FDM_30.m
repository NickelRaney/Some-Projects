%% FDM 30 均匀撒点法主程序 输入N和eps 返回时间、迭代次数、误差、残量
function [tT,cnt,err,res]=main_FDM_30(N,e)
s3=sqrt(3);
tic
%% 构造大的均匀网格
Lout=2;
h=1/N;
x_grid=0:h:Lout;
y_grid=0:h:Lout;
[X,Y]=meshgrid(x_grid,y_grid);
N_tmp=length(x_grid);
%% 标记在正方形内部的点 并生成Ny L R等量
Z=if_in(X,Y);
N_y=N_tmp;
L_x=zeros(N_tmp,1);
R_x=zeros(N_tmp,1);%x轴上无点满足 故从i=2开始
for i=2:N_tmp
    if find(Z(i,:),1)
        L_x(i)=find(Z(i,:),1);
        R_x(i)=find(Z(i,:),1,'last');
    else
        N_y=i-1;
        break;
    end
end
%% 参数及初值
eps=e;
%真解
u_ana=((s3*X+Y-s3/2).*(s3*X+Y-s3/2-2).*(X-s3*Y-1/2).*(X-s3*Y+3/2));
u_ana=u_ana.*Z;

f=ana_sol(X,Y,eps).*Z*h^2;%外力
u=Z*0;%初值
%% 生成每次求解要用到的矩阵和向量
Ai=cell(N_y,1);%每行的矩阵
FLAG=false;%判断最后一行是不是只有一个元素 即是否有落单顶点
%% 处理水平边界
for i=2:N_y
    len=R_x(i)-L_x(i)+1;
     if len==1
         FLAG=true;
         continue;
    end

    dia_x=2*diag(ones(1,len));
    dia_U=-ones(1,len-1);
    dia_L=-ones(1,len-1);
    %% 修正左边的
    x=X(i,L_x(i));
    y=Y(i,L_x(i));
    xL1=(s3/2-y)/s3; xL2=(s3*y-3/2);
    xL=max(xL1,xL2);
    delta=(x-xL)/h;%计算delta
  
    dia_x(1,1)=2/delta;
    dia_U(1)=-2/(delta+1);
    %% 修正右边的
    x=X(i,R_x(i));
    y=Y(i,R_x(i));
    xR1=(2+s3/2-y)/s3; xR2=(1/2+s3*y);
    xR=min(xR1,xR2);
    delta=(xR-x)/h;

    dia_x(end,end)=2/delta;
    dia_L(end)=-2/(delta+1);
    Ai{i}=dia_x+diag(dia_L,-1)+diag(dia_U,1);
end
delta_up=cell(N_y,1);% y方向的情况存在delta_up里  delta 的长度和 u相同 能这么写用到了u在合理的坐标外面是0
delta_down=cell(N_y,1);% 
%% 处理垂直边界
for i=2:N_y
    L=L_x(i);
    R=R_x(i);
    len=R_x(i)-L_x(i)+1;
    if len==1%是否有落单顶点
        continue;%忽略落单顶点
    end
    dia=2*ones(1,len);
    %% 初始化delta_up delta_down
    d_up=Z(i,:).*ones(1,N_tmp);
    d_down=d_up;
    %% 向上修正
    if i==2
        delta_down{i}=d_down;
    else
        %% 向上修正 如果触碰到了上边界 此时不会碰到下边界 因此delat_down的值需要修改
        %% 右边
        if i==N_y||(FLAG&&i==N_y-1)
            st=L;
        else
            st=max(R_x(i+1)+1,L);
        end
        for k=R:-1:st
            [flag,del]=if_up(X(i,k),Y(i,k),h);
            assert(flag||(FLAG&&i==N_y-1));
            d_down(k)=del;
            dia(end-(R-k))=2/del;
        end
       %% 左边
        if i==N_y||(FLAG&&i==N_y-1)
            ed=R;
        else
            ed=min(L_x(i+1)-1,R);
        end
        for k=L:ed
            [flag,del]=if_up(X(i,k),Y(i,k),h);
            assert(flag||(FLAG&&i==N_y-1));
            d_down(k)=del;
            dia((k-L)+1)=2/del;
        end
        delta_down{i}=d_down;
    end
    %% 向下修正 触碰到了下面
    if i==N_tmp
        delta_up{i}=d_up;
    else
       %% 右边
        if i==2
            st=L;
        else
            st=max(R_x(i-1)+1,L);
        end
        for k=R:-1:st
            [flag,del]=if_down(X(i,k),Y(i,k),h);
            assert(flag);
            d_up(k)=del;
            dia(end-(R-k))=2/del;
        end
       %% 左边
        if i==2
            ed=R;
        else
            ed=min(L_x(i-1)-1,R);
        end
        for k=L:ed
            [flag,del]=if_down(X(i,k),Y(i,k),h);
            assert(flag);
            d_up(k)=del;
            dia((k-L)+1)=2/del;
        end
        delta_up{i}=d_up;
    end
    
    Ai{i}=Ai{i}+eps*diag(dia);%加到矩阵上
end

%% 开始迭代!!
Res=0;%先计算初始残量
for i=2:N_y
    if i==N_y
        if R_x(i)==L_x(i)%是否有落单顶点
            continue;%忽略落单顶点
        end
    end
    b=f(i,L_x(i):R_x(i))+eps*(2./(1+delta_up{i}(L_x(i):R_x(i)))).*u(i+1,L_x(i):R_x(i))...
        +eps*(2./(delta_down{i}(L_x(i):R_x(i))+1)).*u(i-1,L_x(i):R_x(i));
    vec=u(i,L_x(i):R_x(i))';
    r=Ai{i}*vec-b';
    Res=Res+norm(r,2)^2;
end
Res=sqrt(Res);

Res0=Res;
cnt=0;
%% Line GS
while Res/Res0>1e-6
    cnt=cnt+1;
    for i=2:N_y
        if i==N_y
            if R_x(i)==L_x(i)
                continue;%忽略落单顶点
            end
        end
        b=f(i,L_x(i):R_x(i))+eps*(2./(1+delta_up{i}(L_x(i):R_x(i)))).*u(i+1,L_x(i):R_x(i))...
            +eps*(2./(delta_down{i}(L_x(i):R_x(i))+1)).*u(i-1,L_x(i):R_x(i));
        vec=Ai{i}\b';
        u(i,L_x(i):R_x(i))=vec';
    end
    Res=0;%计算残差
for i=2:N_y
    if i==N_y
        if R_x(i)==L_x(i)
            continue;
        end
    end
    b=f(i,L_x(i):R_x(i))+eps*(2./(1+delta_up{i}(L_x(i):R_x(i)))).*u(i+1,L_x(i):R_x(i))...
        +eps*(2./(delta_down{i}(L_x(i):R_x(i))+1)).*u(i-1,L_x(i):R_x(i));
    vec=u(i,L_x(i):R_x(i))';
    r=Ai{i}*vec-b';
    Res=Res+norm(r,2)^2;
end
    Res=sqrt(Res);
end
res=Res;
tT=toc;
err=norm(u-u_ana,2)*h;
end
%判断向下有没有出边界 返回delta
function [flag,delta]=if_down(x,y,h)
flag=true;

s3=sqrt(3);
y_d1=(x-1/2)/s3;
y_d2=s3/2-s3*x;
y_d=max(y_d1,y_d2);
if y_d<y-h
    flag=false;
    return;
else
    delta=(y-y_d)/h;
end
end
%判断向上有没有出边界 返回delta
function [flag,delta]=if_up(x,y,h)
flag=true;

s3=sqrt(3);
y_up1=(x+3/2)/s3;
y_up2=(2+s3/2-s3*x);
y_up=min(y_up1,y_up2);
if y_up>y+h
    flag=false;
   % disp('sp')
    delta=(y-y_up-h)/h;
    return;
else
    delta=(y_up-y)/h;
end
end
%判断是否在正方形里面
function rst=if_in(X,Y)
s3=sqrt(3);

rst=(s3*X+Y>s3/2)&(s3*X+Y<s3/2+2)&(X-s3*Y<1/2)&(X-s3*Y>-3/2);
end