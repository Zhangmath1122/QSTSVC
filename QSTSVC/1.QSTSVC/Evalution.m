function [Results] = Evalution(X,Y,para)
% function [Results] = Evalution(X,Y,lamda,C,p)
tic
lamda=para(1);
C=para(2);
p=para(3);
PY = QSTSVC_pre(X,Y,lamda,C,p);%
Result = test(Y,PY);
time=toc;

Results=[Result,time];%