function [x,dx] = satxyz(eph,tcur)
% SYNOPSIS
% 	[position velocity] = satxyz(ephemeris,current_time)
%
% DESCRIPTION
% 	Transform Keplerian elements into inertial 
%	rectangular coordinates.
%
% SEE ALSO
%	satlos
%
% AUTHOR
% 	J.F. Hunzinger 15/7/96

% constants
we = 0.72921151467e-4; % earth rotation rate wrt inertial space WGS-84 [rad/s]
%we = 0; % no rotation (debugging)
UGP = 3.986005e14; % universal graitational potential constant

% -----------------------------------------------------------------------
% YUMA format ephemeris data
% -----------------------------------------------------------------------
no      = eph(1); % 1 - satellite number
status  = eph(2); % 2 - health
e       = eph(3); % 3 - orbit eccentricity
tapl    = eph(4); % 4 - time of applicability
refi    = eph(5); % 5 - inclination angle of orbital plane at reference time
domega  = eph(6); % 6 - rate of right ascension
roota   = eph(7); % 7 - sqrt(semi-major axis of ellipse)  
refomega= eph(8); % 8 - longitude of right ascending node at weekly epoch
w       = eph(9); % 9 - argument of perigree
refM    = eph(10);% 10- mean anomaly at reference time
af0     = eph(11);% 11-
af1     = eph(12);% 12-  
teph    = eph(13);% 13- ephemeris data time (week)

% -----------------------------------------------------------------------
% other directly calculated parameters
% -----------------------------------------------------------------------
a = roota^2; 

% -----------------------------------------------------------------------
% Calculate mean anomally from mean motion
% -----------------------------------------------------------------------
% time since epoch reference
dt = tcur - tapl;
% mean motion, n = sqrt(u/a^3)
n = sqrt(UGP/a^3);
% adjusted mean anomally
M = refM + n*dt;

% -----------------------------------------------------------------------
% Solve Kepler's equation   M = E - e sinE  
% -----------------------------------------------------------------------
% use approximation to E^4
e0 = atan(e*sin(M))/(1-e*cos(M));
e1 = sqrt(1-2*e*cos(M)+e^2);
e2 = e0+asin((-1*e*sin(M+e0))^3/(6*e1));
E = M+e2;


% -----------------------------------------------------------------------
% Correct the inclination
% -----------------------------------------------------------------------
i = refi;

% -----------------------------------------------------------------------
% Correct the right ascension (note we is rotation rate of earth)
% -----------------------------------------------------------------------
omega = refomega + (domega-we)*dt-we*tapl;

% -----------------------------------------------------------------------
% Calculate the radius
% -----------------------------------------------------------------------
r = (1-e*cos(E))*a;

% -----------------------------------------------------------------------
% calculate the true anomaly
% -----------------------------------------------------------------------
denom = 1-e*cos(E);
sinf = sqrt(1-e^2)*sin(E)/denom;
cosf = (cos(E)-e)/denom;
f = atan2(sinf,cosf);
% f = atan(sqrt(1-e^2)*sin(E)/(cos(E)-e));

% -----------------------------------------------------------------------
% efficient calculations
% -----------------------------------------------------------------------
coso = cos(omega);
sino = sin(omega);
cosw = cos(w);
sinw = sin(w);
cosi = cos(i);
sini = sin(i);

Rxq = [
	coso*cosw-sino*cosi*sinw	-coso*sinw-sino*cosi*cosw	sino*sini
	sino*cosw+coso*cosi*sinw	-sino*sinw+coso*cosi*cosw	-coso*sini
	sini*sinw			sini*cosw			cosi
];

q = [r*cos(f) r*sin(f) 0]';
dq = (n*a/sqrt(1-e^2)) * [-sin(f) e+cos(f) 0]';

% -----------------------------------------------------------------------
% transform
% -----------------------------------------------------------------------
x = Rxq*q;
dx = Rxq*dq;
