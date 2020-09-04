%改编自孔鼎问同学之前的代码
%M为刚度矩阵 M = M1 + eps*M2
%为稀疏矩阵形式
function [M,M1,M2]=assembleMatrix_new3(N,eps,w)
    c=cos(w);
    s=sin(w);
    
    alpha1=(2+2*c*s);
    beta1=-(c*s+c^2);
    gamma1=-(s^2+c*s);
    delta1=c*s;
    
    T1=sparse(ones(N-1,1));
    T2=sparse(ones(N-2,1));
    
    A1=diag(T2*beta1,-1)+diag(T2*beta1,1)+alpha1*speye(N-1);
    B1=diag(T2*delta1,1)+gamma1*speye(N-1);
    
    M1=kron(speye(N-1),A1)+kron(diag(T2,1),B1)+kron(diag(T2,-1),B1');
    
    alpha2=(2-2*c*s);
    beta2=-(s^2-c*s);
    gamma2=-(c^2-c*s);
    delta2=-c*s;
    
    T1=sparse(ones(N-1,1));
    T2=sparse(ones(N-2,1));
    
    A2=diag(T2*beta2,-1)+diag(T2*beta2,1)+alpha2*speye(N-1);
    B2=diag(T2*delta2,1)+gamma2*speye(N-1);
    
    M2=kron(speye(N-1),A2)+kron(diag(T2,1),B2)+kron(diag(T2,-1),B2');
    
    M = M1 + eps*M2;
end
    
    
    