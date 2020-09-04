%��ͼ
load('FDM_FINAL_YEAH.mat')
% 
% n=5:8;
% 
% N=2.^n;
% e=[1,1e-2,1e-4,1e-6,1e-8];
% 
% n=length(e);
% m=length(N);
% for i=1:m
%     for j=n:-1:1
%        % [tT(j,i),cnt(j,i),err(j,i),res(j,i)]=main_FDM_60(N(i),e(j));
%         str=sprintf('N=%d, eps=%g. cnt=%d err=%g. using %g s. res=%g',N(i),e(j),cnt(j,i),err(j,i),tT(j,i),res(j,i));
%         disp(str)
%     end
% end


%% ��ͬeps ����N
figure(1)
for i=1:n
    str=sprintf('eps=%g',e(i));
    plot(log2(N),log2(err(i,:)),'Displayname',str);
    hold on
end
xlabel('log_2(N)')
ylabel('log_2(err)')
%% ��ͬeps ����������N
figure(2)
for i=1:n
    str=sprintf('eps=%g',e(i));
    plot(log2(N),log2(cnt(i,:)),'Displayname',str);
    hold on
end
xlabel('(N)')
ylabel('iteration time')

%% ��ͬN ����������eps
figure(3)
for i=1:m
    str=sprintf('N=%d',N(i));
    plot(log10(e),log10(cnt(:,i)),'Displayname',str);
    hold on
end
xlabel('log_{10}(eps)')
ylabel('log_{10}(iteration time)')