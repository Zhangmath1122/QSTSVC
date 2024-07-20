function [PY] = QSTSVC_pre(X,Y,lamda,C,p)
%% 数据X,Yreal,lamda,C,p
[m,n]=size(X);
Y0=NNG(X,Y,p);     %初始标签
K=max(Y);%分类数K
totalZ=zeros(K,(n^2+3*n+2)/2);%最终Z矩阵（K个超曲面）
%% 数据处理
S=zeros(m,(n^2+3*n+2)/2); %数据拉伸后的矩阵S
for i=1:m  %m个数据
    S1=[];
    for a=1:n
        for b=a:n
            if a==b
                s=0.5;
            else
                s=1;
            end
            s1=s*X(i,a)*X(i,b);
            S1=[S1 s1];
        end
    end
    S1data{i}=S1;
    S(i,:)=[S1data{i} X(i,:) 1];
end
prev = 0;
prevY = 0;
diff = 1;
diff2 = 1;
times = 0;
PY=Y0;
while  ((norm(diff) > 0.1) && (norm(diff2) ~= 0) && times <= 5)
    times = times + 1;

    %数据处理
    for i=1:K
        M1=S(PY==i,:);  %第i类数据矩阵M1
        M2=S(PY~=i,:);  %其他类数据矩阵M2
        if (isempty(M1))
            continue;
        end
        if (isempty(M2))
            continue;
        end
        A1=X(PY==i,:); %第i类原始数据矩阵A
        Z0=qslssvr_pre(A1(:,1:n-1),A1(:,n),C);
%         Za=FirstStep(A1);
%         Z0=[zeros((n^2+n)/2,1);Za];
%         Z0=FirstStep(M1);
%         Z0=Z0(1:(n^2+3*n+2)/2);
        m1=size(M1,1);  %第i类数据个数m1
        m2=size(M2,1);  %其他类数据个数m2
        G=zeros((n^2+3*n+2)/2,(n^2+3*n+2)/2); %构造G
        for j=1:m1
            G=G+M1(j,:)'*M1(j,:);  %得到G
        end
        %%求解
        tic
        ite=0;
        som=1;
        tol=0.000001;
        while som>tol && ite<10
            ite=ite+1;
            Z=Z0;
            o1=zeros((n^2+3*n+2)/2,m2);
            o2=zeros(m2,m2);
            H=[G o1;o1' o2];%构造G
            H=H+1e-6*eye(size(H,1)); %加一个很小的单位阵
            H=0.5*(H+H');%对角化处理
            f=[zeros((n^2+3*n+2)/2,1);lamda*ones(m2,1)];
            A=[-diag(sign(M2*Z))*M2,-eye(m2)];
            b1=-1*ones(m2,1);
            lb=[-inf*ones((n^2+3*n+2)/2,1);zeros(m2,1)];
            gamma=quadprog(H,f,A,b1,[],[],lb,[],[],optimset('display','off'));
            Z0=gamma(1:(n^2+3*n+2)/2);
            som=norm(Z-Z0);
        end
        totalZ(i,:)=Z0';
    end
    [d, pY] = min(abs(S * totalZ'), [], 2);
    PY = pY;
    diff = sum(d) - prev;
    prev = sum(d);
    diff2 = pY - prevY;
    prevY = pY;
end



