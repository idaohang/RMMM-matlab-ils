function [wgs] = xyz2wgs(pos)
% SYNOPSIS
% 	[wgs84] = xyz2wgs(X)
%
% DESCRIPTION
%       Convert the ECEF X-Y-Z Coordinates to WGS-84 
%	Latitude, Longitude and Height Coordinates.
%       Here X = [x y z] and wgs84 = [lat lon h].
%
% SEE ALSO
%	dxyz2ned, satxyz	
%
% AUTHOR
% 	J.F. Hunzinger 17/7/96

esma = 6378137; % earth semi-major axis
ee2 = 0.00669437999013; % earth eccentricity squared

x = pos(1);
y = pos(2);
z = pos(3);

p = 1e-30;
% protect against zero
if(abs(x) < p),
  x = p * sign(x); 
end;
if(abs(y) < p),
  y = p * sign(y);
end;
if(abs(z) < p),
  z = p * sign(z);
end;

denom = sqrt(x^2+y^2);
long = atan2(y,x);
lat = atan(z/((1-ee2)*denom));
% iterate here for better precision 
% see CMC XYZGEO()

% radius of curvature along meridian line
rp = esma/sqrt(1-ee2*sin(lat)*sin(long));
clat = cos(lat);
if(abs(clat) < 1e-14),
  height = abs(z)-(1-ee2)*rp;
else
  height = denom/clat-rp;
end;

wgs = [lat long height]';