%% 按照题目要求的自动求解程序
clear; clc;
w = pi/4;
c = cos(w); s = sin(w);   
F = @(x,y)(-2*(c^4*y^2+s^3*y*(-1+s*y)+c*s*(1+3*s*x)*(-1+2*s*y)-... 
   c^3*(y+6*s*x*y)+c^2*s*(3*x+6*s*x^2+2*y*(1-2*s*y))))+eps*(...
   -2*(c^4*x^2+s^3*x*(1+s*x)-c*s*(1+2*s*x)*(-1+3*s*y)+... 
   c^3*x*(-1+6*s*y)-c^2*s*(2*x+4*s*x^2+3*y*(1-2*s*y))));
for N = [32,64,128,256]
    for eps = [1,1e-2,1e-4,1e-6,1e-8]
        [u,k,~,iniTime,solTime] = FDM45solver(N,eps,F);

        x = -1:1/N:1; 
        y = 0:1/N:2;
        [X,Y] = meshgrid(x,y); 
        X = X/sqrt(2); Y = Y/sqrt(2);
        v = (Y*cos(w)-X*sin(w)).*(-Y*cos(w)+X*sin(w)+1).*...
            (X*cos(w)+Y*sin(w)).*(-X*cos(w)-Y*sin(w)+1);    % 真解
        T = (Y*cos(w)-X*sin(w)>=0)&(-Y*cos(w)+X*sin(w)+1>=0)&...
            (X*cos(w)+Y*sin(w)>=0)&(-X*cos(w)-Y*sin(w)+1>=0);
        v = T.*v;
        
%         mesh(X,Y,u)
        err = norm(u(:)-v(:),2)/N;

        fprintf('w=pi/4\tN=%d\teps=%e\terr=%e\titers=%d\tsolTime=%e\tiniTime=%e\n',...
            N,eps,err,k,solTime,iniTime);

    end
end