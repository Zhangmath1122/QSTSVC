function [Z0]=qslssvr_pre(X,Y,C)
lamda=C;
X=X';
[m,n]=size(X); %训练集input的维数
W=rand(m,m); %初始化W矩阵
W1=triu(W);  %取上三角矩阵
[l,h]=find(W1');
Q=[l';h';1:((m^2+m)/2)];
M=zeros(m,(m^2+m)/2);
for i=1:1:n
    for k=1:m
        for j=1:m
            for p=1:(m^2+m)/2
                if (Q(1,p)==j&&Q(2,p)==k)||(Q(1,p)==k&&Q(2,p)==j)
                    M(j,p)=X(k,i);
                    M(k,p)=X(j,i);
                end
            end
        end
    end
    data{i}=M'; %构造M矩阵
end
U=zeros((m^2+3*m+2)/2,n);
for k=1:n
    A_x=[];
    for i=1:m
        for j=i:m
            if i==j
                s=0.5;
            else
                s=1;
            end
            t=s*X(i,k)*X(j,k);
            A_x=[A_x t];
        end
    end
    AList=A_x;
    L=X(:,k);
    Ldata{k}=L';
    ALdata{k}=AList;
    U(:,k)=[ALdata{k} Ldata{k} 1]';
end
Q=lamda*U*U';
r_term=lamda*U*Y;
H_nonlinear=zeros((m^2+m)/2,(m^2+m)/2);
for i=1:n
    H_nonlinear=H_nonlinear+data{i}*data{i}'; %计算矩阵H_nonlinear
end
H_linear=diag(2*ones(m,1));
H=blkdiag(H_nonlinear,H_linear,0);
ll_term=H+Q;
I=eye((m^2+3*m+2)/2);
l_term=ll_term+1e-6*I;
z=l_term\r_term;

v=z(1:(m^2+m)/2);
D=z((m^2+m+2)/2:(m^2+3*m)/2);
c=z((m^2+3*m+2)/2);
for i=1:m
    for j=i:m
        index=sum(m:-1:m-i+2)+j-i+1;
        ww(i,j)=v(index);
    end
end
W=ww'+triu(ww,1);
D1=D;
c1=c;
o1=zeros(m,1);
W2=[W o1;o1' 0];
D2=[D1;-1];
D2=D2';
W3=triu(W2);  %取上三角矩阵
W4=[];
for i=1:m+1
    for j=i:m+1
        t1=W3(i,j);
        W4=[W4 t1];
    end
end
Z0=[W4,D2,c1];
Z0=Z0';