function [e] = elevangl(rec,sat)
% SYNOPSIS
% 	[eauv] = elevangl(receiver_ecef_xyz,sat_ecef_xyz)
%
% DESCRIPTION
%	Calculates the unit vector in the direction of 
%       the satellite from the receiver. The coordinate
%       system used is ECEF XYZ.
%
%	Returns: eauv = Elevation angle unit vector
%
% SEE ALSO
%
% AUTHOR
% 	J.F. Hunzinger 21/6/96

e = sat - rec;
% normalize
e = e / sqrt(sum(e.^2));
