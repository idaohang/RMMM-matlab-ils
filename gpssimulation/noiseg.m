function [d] = noiseg(m,v)
% SYNOPSIS
% 	noiseg(mean,variance)
%
% DESCRIPTION
%	Return normal random number with mean m, variance v.
%
% AUTHOR
% 	J.F. Hunzinger 21/6/96

d = randn(1)*v+m;
%d=0;