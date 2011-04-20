function [a] = phasplit(phase)
% SYNOPSIS
% 	[parts] = phasplit(phase)
%
% DESCRIPTION
%	Calculate the fractional and integer parts of
%       carrier phase (not range = phase*wavelength)
%
%	parts = [fractional_part integral_part]
%
% SEE ALSO
%
% AUTHOR
% 	J.F. Hunzinger 21/6/96

%if phase<0,
%  n = ceil(phase);
%else
  n = floor(phase);
%end;
f = phase - n;
a = [f n];