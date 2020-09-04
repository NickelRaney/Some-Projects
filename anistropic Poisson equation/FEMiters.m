clear; clc;
N = 64;

iterVec = zeros(9,2);
for w = [0,pi/6,pi/4]
    c = cos(w); s = sin(w);
    F = @(x,y)(-2*(c^4*y^2+s^3*y*(-1+s*y)+c*s*(1+3*s*x)*(-1+2*s*y)-... 
       c^3*(y+6*s*x*y)+c^2*s*(3*x+6*s*x^2+2*y*(1-2*s*y))))+eps*(...
       -2*(c^4*x^2+s^3*x*(1+s*x)-c*s*(1+2*s*x)*(-1+3*s*y)+... 
       c^3*x*(-1+6*s*y)-c^2*s*(2*x+4*s*x^2+3*y*(1-2*s*y))));
    i = 0;
    for eps = [1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8]
        i = i+1;
        [~,k,~,~,~] = FEMsolver(N,eps,w,F);
        iterVec(i,:) = [-log10(eps),k];
    end
    plot(iterVec(:,1),iterVec(:,2));
    hold on
end
grid on
title('N=64时迭代次数与\epsilon关系')
xlabel('-log10(\epsilon)')
ylabel('迭代次数')
legend('\omega=0','\omega=\pi/6','\omega=\pi/4')