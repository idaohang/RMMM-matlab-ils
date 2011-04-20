function [ea,azmth] = satlos(xr,xs)
% SYNOPSIS
% 	[elevangl azmth] = satlos(recv_pos_ecef,sat_pos_ecef);
%
% DESCRIPTION
%	Calculate the elevation angle of the satellite from
%	the receiver point of view and also the azimuth angle.
%
% SEE ALSO
%	satxyz
%
% AUTHOR
% 	J.F. Hunzinger 17/7/96

% vector to satellite
dx = xs-xr;

% receiver lat and long (and height)
wgs84 = xyz2wgs(xr);

% calculate N,E,DOWN of the vector to the satellite
ned = dxyz2ned(dx,wgs84);

p = 1e-30;
% protect against zero
for i=1:3,
  if(abs(ned(i)) < p),
    ned(i) = p * sign(ned(i)); 
  end;
end;

n = ned(1);
e = ned(2);
d = ned(3);

% from the receiver's viewpoint
% find the elevation angle (may be used to see if the line-of-site
% is below the horizon.
ea = atan(-d/sqrt(n^2+e^2));

% find the azimuth angle (may be used to determine if an object
% is blocking the line-of-site.
azmth = atan2(e,n);

