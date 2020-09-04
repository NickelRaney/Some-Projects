%% 按照题目要求的自动求解程序
clear; clc;

for w = [0,pi/6,pi/4]
    c = cos(w); s = sin(w);   
F = @(x,y)(-2*(c^4*y.^2+s^3*y.*(-1+s*y)+c*s*(1+3*s*x).*(-1+2*s*y)-... 
   c^3*(y+6*s*x.*y)+c^2*s*(3*x+6*s*x.^2+2*y.*(1-2*s*y))))+eps*(...
   -2*(c^4*x.^2+s^3*x.*(1+s*x)-c*s*(1+2*s*x).*(-1+3*s*y)+... 
   c^3*x.*(-1+6*s*y)-c^2*s*(2*x+4*s*x.^2+3*y.*(1-2*s*y))));
    for N = [32,64,128,256]
        for eps = [1,1e-2,1e-4,1e-6,1e-8]
            %             [u,k,~,iniTime,solTime] = FEMsolver(N,eps,w,F);
            [u,k,~,iniTime,solTime] = FDMsolver(N,eps,w,F);
            
            x = 0:1/N:1; y = 0:1/N:1;
            [X,Y] = meshgrid(x,y);
            X1 = cos(w)*X-sin(w)*Y;
            Y1 = sin(w)*X+cos(w)*Y;

            v = (Y1*cos(w)-X1*sin(w)).*(-Y1*cos(w)+X1*sin(w)+1).*...
                (X1*cos(w)+Y1*sin(w)).*(-X1*cos(w)-Y1*sin(w)+1);    % 真解
            err = norm(u(:)-v(:),2)/N;
            
            fprintf('w=%e\tN=%d\teps=%e\terr=%e\titers=%d\tsolTime=%e\tiniTime=%e\n',...
                w,N,eps,err,k,solTime,iniTime);
            
%             %% 绘制图像
%             % 残差图像
%             figure(1);
%             semilogy(resVec);
%     
%             % 解的图像
%             figure(2);
%             x = 0:1/N:1; y = 0:1/N:1;
%             [X,Y] = meshgrid(x,y);
%             X1 = cos(w)*X-sin(w)*Y;
%             Y1 = sin(w)*X+cos(w)*Y;
%             mesh(X1,Y1,u);
%     
%             % 误差的图像
%             figure(3)
%             mesh(X1,Y1,u-(Y1*cos(w)-X1*sin(w)).*(-Y1*cos(w)+X1*sin(w)+1).*...
%                 (X1*cos(w)+Y1*sin(w)).*(-X1*cos(w)-Y1*sin(w)+1));
    
        end
    end
end