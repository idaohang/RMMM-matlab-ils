function [cl] = rclock(cl,st,R)
% DESCRIPTION: Simulation of receiver clock bias
% AUTHOR:      Valeri Perepetchai 03/00.

A=[1 st;0 1];
cl=A*cl+R*randn(2,1); % in sec
%cl(1)=abs(cl(1));
cl=abs(cl);
