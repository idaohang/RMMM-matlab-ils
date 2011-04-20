function [Rup] = rcfac(st)
% DESCRIPTION: Computing Chol factor of cov. matrix of clock bias
% AUTHOR:      Valeri Perepetchai 03/00.

h0=9.4e-20;
h1=1.8e-19;
h2=3.8e-21;
c11=h0/(2*st) + 2*h1*st^2+(2/3)*pi^2*h2*st^3;
c12=2*h1*st+pi^2*h2*st^2;
c22=h0/(2*st) + 2*h1+(8/3)*pi^2*h2*st;
Rup=chol([c11 c12;c12 c22]);% upper triang. Chol. factor
