function [pos] = rcvrmdl(h,r,v,epoch);
% SYNOPSIS
% 	[pos] = rcvrmdl(height,radius,epoch)
%
% DESCRIPTION
%	Calculate the position of rover receiver
%	at each epoch
%
%       height = height of aircraft above the landing site
%                the landing site is assumed to be (EARTH_RADIUS 0 0)
%       radius = radius of the aircraft circular flight path
%       epoch = time instant of flight
%      
%       
%       all measurements are in ECEF-XYZ coordinates 
%
% AUTHOR
% 	J.F. Hunzinger 13/11/96
%

EARTH_RADIUS 	= 6378137; 	% earth radius

% calculate the angular velocity
a = v/r;

% calculate the current angular position
p = rem(a*epoch,360);

% calculate the position
x = EARTH_RADIUS+h;
y = r*sin(p);
z = r*cos(p);
pos = [x y z];

% calculate the velocity
% not required currently

