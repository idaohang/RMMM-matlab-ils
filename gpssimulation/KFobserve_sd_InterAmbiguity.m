function [visible,sd,iNs] = KFobserve_sd_InterAmbiguity(x_rcvr1,x_rcvr2,dx_rcvr2e,sigma,t,min_elev,ephd,n_vis_sat_max)
% DESCRIPTION
%       Determine the visible satellites and the range from each receiver
%       to the satellites (also calculates the difference in those ranges
%       for each satellite for efficiency).
%

L1_F 		= 1575.42e6; 	% L1 carrier frequency
C 		= 3e+8; 	% speed of light
L1_W 		= C/L1_F; 	% L1 wavelength


% determine the estimated mid baseline point
mid_base_line = ((dx_rcvr2e)/2 + x_rcvr1)';

[n_sat n_par] = size(ephd);

% calculate the satellite positions and elevation angles
for i=1:n_sat,
  % calculate satellite position
  sv(i,:) = satxyz(ephd(i,:),t)';
  % check if visible
  [ea,az] = satlos(mid_base_line,sv(i,:)');
  vis(i) = ea;
end;

% select the visible satellites
[viss,index] = sort(vis); % increasing order
for i=1:n_sat,
  if(viss(i) > min_elev), break; end;
end;
visible = index(n_sat:-1:i);
n_vis_sat = length(visible);
% sort the non-reference satellites according to satellite number (see below)
visible = [visible(1) visible(2) sort(visible(3:n_vis_sat))];

if(n_vis_sat > n_vis_sat_max),
  visible = visible(1:n_vis_sat_max);
  n_vis_sat = n_vis_sat_max;
end;
%ref = index(n_sat);

% compute the actual noisy ranges and delta ranges, sd
svv = sv(visible,:);
for i=1:n_vis_sat,
  % compute the ranges
  r_rcvr1(i) = sqrt(sum((svv(i,:)-x_rcvr1).^2));
 
  r_rcvr2(i) = sqrt(sum((svv(i,:)-x_rcvr2).^2));
  

  tmp = phasplit((r_rcvr1(i)-r_rcvr2(i))/L1_W);
  sd(i) = tmp(1)+noiseg(0,sigma)/L1_W;
  iNs(i)= tmp(2);  
end;












