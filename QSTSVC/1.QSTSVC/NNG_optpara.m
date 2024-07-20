function [optp] = NNG_optpara(X,Y)
p=(1:5);
Scope={p};%
ParaNum=length(Scope);             %需要调参的个数   Scope元胞数组，每个元包对应一个参数的网格向量
TestNum=1;
TestEach=zeros(ParaNum,1);         %循环完返回每个参数网格的个数
for ii=1:ParaNum
    TestEach(ii)=length(Scope{ii});  %每个参数的网格个数
    TestNum=TestNum*TestEach(ii);    %计算所有参数组合个数
end

% Turn training mesh into a table
[P,S]=deal(ones(1,ParaNum),cell(1,ParaNum));  %表格初始大小   S:元包放每组参数  P:向量放对应的参数组合的网格位置
T=cell(TestNum,1);                            %返回所有的参数组合存入元包，每个元包存一组参数值
for ii=1:TestNum                              %循环第ii组参数组合
    TN=ii-1;                                      %临时变量
    for pp=1:ParaNum
        P(pp)=rem(TN,TestEach(pp))+1;             %第ii个参数组合中第pp个参数的位置
        tmp=Scope{pp};                            %第pp个参数网格向量
        S{pp}=tmp(P(pp));                         %给S元包赋每组参数值
        TN=floor(TN/TestEach(pp));                %构造的临时变量用来计算参数组合的网格位置
    end
    T{ii}=S;                                      %将第ii组参数组合赋给T{ii}
end
R=T;                                        %返回所有的参数组合构成的元包
s=length(R);
Results1=zeros(s,1);
%计算系数
for j=1:s
    p = R{j}{1,1};
    [PY] = NNG(X,Y,p);
    Result = test(PY,Y);
    Results1(j,1)=Result(2);
    disp(['当前迭代次数：' num2str(j)])
end
[~,b]=max(Results1);
optp=R{b}{1,1};
% para=[optlamda optC optp];