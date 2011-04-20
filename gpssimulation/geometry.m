function [E] = geometry(x_rcvr1,dx_rcvr2e,visible,ephd,t)
% SYNOPSIS
% 	[E] = 
%          geometry(x_rcvr1,dx_rcvr2e,visible,ephd,time)
%
% DESCRIPTION
%       Determine the geometry matrix given an estimated baseline.
%
% SEE ALSO
%
% AUTHOR
% 	Jason F. Hunzinger 17/12/96
%
% NOTES:
%	- all distances are measured in meters
%	- all times are measured in seconds
%

% calculates:
%	n_vis_sat	number of visible satellites = length(visible)
%       E               elevation unit vectors


% determine the estimated mid baseline point
mid_base_line = ((dx_rcvr2e)/2 + x_rcvr1)';

n_vis_sat = length(visible);

% calculate the satellite positions and elevation angles
for i=1:n_vis_sat,
  % calculate satellite position
  sv(i,:) = satxyz(ephd(visible(i),:),t)';
  E(i,:) = elevangl(mid_base_line',sv(i,:));
end;











