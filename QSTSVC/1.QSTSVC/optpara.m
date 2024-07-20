function [Results1,para]= optpara(X,Y)

[optp] = NNG_optpara(X,Y);
[Results1,para] = QSTSVC_optpara(X,Y,optp);