% function [ecco,inco,w_po,capomo,L_po,obl,azo,wo,fo,Mo,ao,t,obl_so] = orbelm_swift(t);
function [ecco,inco,w_po,capomo,L_po,obl,wo,fo,Mo,ao,Lso,t] = orbelm_swift(t)
%% HEADER: orbelm_swift.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%Description: Returns orbital elements for Saturn system by interpolation
%              into elements generated by solar system mass / gravity 
%              integration from current epoch to 1e6 years in past.
%
%%Inputs:
%       t: vector of yrs since current epoch 
%          (valid range: -1e6 <= t <= 0)

% MISHA - looks like it has to be row vector, so let check it and convert
% if this is column
if iscolumn(t) 
    t=t.'; 
end

%
%%Outputs:
%       ecco: eccentricty
%       inco: incidence angle [radians]
%       w_po: longitude of perihelion [radians]
%       capomo: longitude of ascending node [radians]
%       L_po: Ls of perihelion [radians]
%       obl: Obliquity of Saturn [radians] (calculated from secular theory)
%       wo: argument of periapse [radians]
%       fo: true anomaly [radians]
%       Mo: mean anomaly [radians]
%       ao: semi-major axis (AU)
%        t: input time vector
%
%%Calls: None
%
%%Required Files: orbelm_swift2.mat
%
%%Source:
%      Integration results of swift orbital position integration code
%
%%Author: Alexander G. Hayes (hayes@gps.caltech.edu)
%
%%Change Log: 
%            09/13/08: Draft Code Escalated to V1
%            01/10/09: Modified Obliquity Calculations and updated Secular
%                      Constants (hayes@gps.caltech.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% BEGIN CODE: orbelm_swift.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%inputs
%t = vector of time in years from J2000
%%outputs (ALL SATURN EPHEMERIS VALUES!)
%ecco = eccentricity
%inco = orbital inclination
%w_po = longitude of perihelion
%capomo = longitude of ascending node
%L_po = Ls of perihelion
%obl = obliquity
%wo = argument of periapse
%fo = true anomaly
%Mo = mean anomaly
%ao = semi-major axis
%t = input time vector
%obl_so = obliquity for constant spin pole
%

%start by giving zero length vectors to outputs (so we don't get warnings)
ecco = [];
inco = [];
w_po = [];
capomo = [];
wo   = [];
fo   = [];
Mo   = [];
ao   = [];
L_po =[];
obl = [];

%in this case we want time to be the negative of input array
t = -t;

%load that data
%matfile = 'orbelm_swift2.mat';
matfile = 'orbelm_swift2_Misha.mat';
if exist(matfile) ~= 2;
    disp(['Can Not Find Required ' matfile ...
        ' Please Ensure It Is In You Matlab Path...']);
    return;
end;
load(matfile,'ecc','inc','w_p','capom','w','f','M','Lm','Ls','a','t0');
names = {'ecc','inc','w_p','capom','w','f','M','Lm','Ls','a','t0'};
num_plan = size(ecc,2);


%make sure int is within the range of t0
ind = find( (t <= max(t0)) & (t >= min(t0)) );

%check to make sure we are in the range of the data
if length(ind) == 0;
    disp(['Problem: Requested Times Not Within Time Period of ' matfile ...
          ' ; t Must Be In Range ' num2str(min(-t0)) ...
          ' <= t <= ' num2str(max(-t0)) ' Returning...']);
    return;
elseif length(ind) < numel(t);
    t = t(ind);
    disp(['Warning: Requested Times Not ALL Within Time Period of ' matfile ...
          ' ; t Must Be In Range ' num2str(min(-t0)) ...
          ' <= t <= ' num2str(max(-t0)) ' Output Truncated...']);
end;

%create the interpolated arrays
num_t = numel(t);
interp_type = 'phcip';
for i=1:length(names);
    %eval([names{i} '= zeros(num_plan,num_t);']);
    for j=1:num_plan;
        eval([names{i} 'o(j,:) = interp1(t0,' names{i} '(:,j),t,interp_type);']);
    end;
end;
wo = mod(wo,2.*pi);
w_po = mod(w_po,2.*pi);


%we want to reverse time again for the secular obliquity varaitions
t = -t;

%calculate the obliquity (WARD 2004, HAMILTON, 2004, WARD 1974, BRETAGNON 1973)
obl_J2000 = 26.73;
%define required constants
n =  [];%Heliocentric mean motion 
w=  [] ;%Spin Frequency of Saturn
J2 = 1.6297e-2; %coefficient of quadrapole moment of gravity field
q = 0.05164; %(French 1993)
    %Effective Quadrupole Moment of satellite sytem 
    %q/J2 = ratio of soalr torque on the satellites to that directly
    %exterted on planet
l = 0.00278; %(French 1993)
lambda = 0.2199; %(Hubbard and Marley (1989)
                 %Moment of Inertia of Saturn normalized to MsR^2
%define precessional constant (WARD 2004 Eq (1))
%alpha = 3/2.*n.^2./w.*(J2+q)./(lambda+l);
alpha = 0.8306/3600; %deg/year
%projected alpha onto mean orbit plane
alpha_cosobj = alpha.*cos(obl_J2000.*pi./180);%0.7380./3600;%0.7427./3600; %deg/year
%large amplitude terms for saturns orbit (WARD 2004, TABLE 1 from BRETAGNON 1973)
gj = [-26.34,-2.99,-0.692]./3600; %deg/yr
Ij = [0.910,0.045,0.064]; %deg
%define phase offset (BRETAGNON 1973)
delta_j = [125.383423,316.173594,201.171559];%deg
%define initial pole azimuth
az_o = 59.9392-90; %59.120 from Ward, 2004, 0.8192 from JPL SPICE analysis
                               %Projection of Spin Vector onto orbit plane
                               %Clockwise angle from vector perpendicular
                               %   to line of ascending node
                               % Instantaneous Value at J2000
%define initial longitude of ascending node
omega_o = 113.6669;
%obl (WARD 2004 EQ.(5) / WARD 1974 EQ. (32)
for i=1:3;
    obl_secular(i,:) = ( gj(i).*Ij(i)./(alpha_cosobj+gj(i)).*...
    (sin( (alpha_cosobj.*t+gj(i).*t+delta_j(i) - az_o - omega_o).*pi./180) - ...
    (sin( (delta_j(i) - az_o - omega_o).*pi./180) ) ) );
end;
obl = obl_J2000 - sum(obl_secular,1);
obl = obl.*pi./180;

%calculate the azimuth change in Titan's Spine Pole (WARD 1974)
%%First WARD 1974 EQ 34
for i=1:3;
    obl_secular2(i,:) = ( gj(i).*Ij(i)./(alpha_cosobj+gj(i)).*...    
    (sin( (delta_j(i) - az_o - omega_o).*pi./180) ) );
end;
obl_o_star = obl_J2000 + sum(obl_secular2,1);
%now define spin axis axzimuth (Ward 1974 EQ 36 and 37)
obl_o = obl_J2000.*pi./180;
for i=1:3;
    obl_secular3(i,:) = ( gj(i).*Ij(i)./(alpha_cosobj+gj(i))./sin(obl_o).*...
    (gj(i).*cos(obl_o) + alpha)./(gj(i)+alpha_cosobj).*...
    (cos( (alpha_cosobj.*t+gj(i).*t+delta_j(i) - az_o - omega_o).*pi./180) - ...
    (cos( (delta_j(i) - az_o - omega_o).*pi./180) ) ) );
end;
pole_lon = az_o - alpha.*t.*cos(obl_o_star.*pi./180) - ...
    capomo.*180./pi + omega_o + sum(obl_secular3,1);
pole_lon = mod(pole_lon.*pi./180,2.*pi);
azo = pole_lon;

%calculate Ls of perihelion (Lsp) using the position of Titans spin pole
L_po = mod(wo-azo-pi/2,2.*pi);

%Obliquity for constant spin pole
obl_o = 26.7.*pi./180; %Obliquity at t=0
inco_o = 2.4848.*pi./180; %Inclination at t=0
az_o = (59.9392-90).*pi./180; %pole azimuth at t=0
omega_o = 113.6669.*pi./180; %Longitude of Ascending Node at t=0
%obl (WARD 1974 EQ (14)
obl_so = acos( ...
    cos(inco).*(cos(inco_o).*cos(obl_o) + sin( inco_o).*sin(obl_o).*sin(az_o) )+ ...
    sin(inco).*sin(capomo-omega_o).*sin(obl_o).*cos(az_o) + ...
    sin(inco).*cos(capomo-omega_o).*(sin(inco_o).*cos(obl_o)-cos(inco_o).*sin(obl_o).*sin(az_o)));
%
%This is my guess for Lso - Misha
%
coeff=57.3;
Lso=Lso/coeff;
%% END CODE: orbelm_swift.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%