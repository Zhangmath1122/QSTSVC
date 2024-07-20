function [PY] = NNG(X,Y,p)

[m,~]=size(X); %数据设定m*n
k=max(Y);%分类数K
D = pdist2(X,X);% 计算样本之间的距离矩阵
% 初始化索引矩阵
I1 = zeros(m,p); %p近邻
I2 = zeros(m,1); %对应距离
% 对于每个点，找到距离最小的p个点
for i = 1:m
    [y2,ind] = sort(D(i,:)); % 按距离升序排序
    I1(i,1:p) = ind(:,2:p+1); % 邻接矩阵
    I2(i,1:p) = y2(:,2:p+1);  % 距离
end
G = graph(I1,repmat((1:m)',1,p));% 创建初始无向图 repmat((1:m)',1,p)
%plot(G,'XData',X(:,1),'YData',X(:,2)); % 画图
L = conncomp (G);% 初始化聚类标签
t = max (L);%无向图目前类别数
while t ~= k
    % 如果当前聚类数目不等于目标聚类数目k，就继续循环
    if t < k % 分裂聚类
        % 如果当前聚类数目小于k，就分裂聚类
        [i1, j1]=find(I2==max(max(I2)));
        if length(i1)>1
            i1=i1(1);j1=j1(1);
        end
        % 找到最大的边权重和对应的两个样本
        j3=I1(i1,j1);% 根据索引找到列号
        % 断开这两个样本之间的连接
        G = rmedge(G,i1,j3);
        %plot(G,'XData',X(:,1),'YData',X(:,2),'ZData',X(:,3)); % 画图
        I2(i1,j1) = 0; %更新距离矩阵
        L = conncomp (G);% 当前聚类标签
    else % 合并聚类
        % 如果当前聚类数目大于k，就合并聚类
        H = zeros(t,t);
        for i=1:t-1
            for j=i+1:t
                dij = D(L==i,L==j);
                H(i,j) = max(max(min(dij,[],1)),max(min(dij,[],2)));
                % Hausdorff distance max(行最小值中的最大值；列最小值的最大值)
                %行最小值min(dij,[],2)
                %列最小值min(dij,[],1)
            end
        end
        H=H+1000000000*tril(ones(t));
        [i2,j2]=find(H==min(min(H)));
        if length(i2)>1
            i2=i2(1);j2=j2(1);
        end
        % 合并这两个聚类，使得原来属于不同的聚类的样本现在属于同一个聚类
        L(L==j2)=i2;
        L(L>j2) = L(L>j2) - 1;
        %t = numel(unique(L));%无向图目前类别
    end
    t = numel(unique(L));%无向图目前类别数
end
PY=L;