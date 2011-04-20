function [ned] = dxyz2ned(dx,wgs)
% SYNOPSIS
% 	[ned] = xyz2ned(dx,wgs84_pos)
%
% DESCRIPTION
%	Convert an ECEF X-Y-Z difference vector
%	to a [North East Down] coordinate vector.
%	To simplify the conversion, it is necessary
%	to provide the wgs84 coordinates of the origin
%	of the difference vector.
%
% SEE ALSO
%	xyz2wgs, satlos
%
% AUTHOR
% 	J.F. Hunzinger 17/7/96


lat = wgs(1);
long = wgs(2);

% compute the xyz to ned transformation matrix
slat = sin(lat);
clat = cos(lat);
slong = sin(long);
clong = cos(long);
A = [
	-clong*slat	-slong*slat	clat
	-slong		clong		0.0
	-clong*clat	-slong*clat	-slat
  ];

ned = A*dx;
