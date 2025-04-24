% function Calcolo_rad_sol_alb_IR_12x24_f_multi_alb_b
%Acc_XYZ_sole
%  Code to calculate the radiation pressure on a Galileo satellite du to direct solar radiation
%  Earth solar albedo and earth IR
%
%  Feb, 23th 2021
%
% -The program start loading the parameters, general and referred to the satellite
% -The time span and the sample rate for simulation is fixed
% -The keplearian parameters are calculated within the simulation period 
%
%
fid1=fopen('Acc_gauss_SRP.out','w');
fid2=fopen('Acc_gauss_AL.out','w');
fid3=fopen('Acc_gauss_IR.out','w');

fid4=fopen('Acc_DYB_DB.out','w');
fid5=fopen('Acc_XYZ_sole.out','w');
fid6=fopen('Acc_XYZ_albedo.out','w');

% carica parametri assenti in David
addpath('C:\Users\Carlo\OneDrive - INAF - Istituto Nazionale di Astrofisica\G4S\Box_wing\Confronto Box Wing Surf Surf_Visco\Script di Confronto\Massimo_risultati\Celestial_computing')
addpath('C:\Users\Carlo\OneDrive - INAF - Istituto Nazionale di Astrofisica\G4S\Box_wing\Confronto Box Wing Surf Surf_Visco\Script di Confronto\Massimo_risultati\');
const=read_par('constant_SI.par');              % carica le costanti SI
astrody=read_par('astrody_SI.par');             % carica parametri astrofisici    
par_bw=read_box_wing_par('box_wing_meta.par');  % carica i parametri del boxwing

astrody.GM = 398600.4418;
%%


%
% Select the satellite on the basis of the initial conditions from ESA
% Galileo Meta data or from the precise orbits provided in the sp3 files
%
% satellite='E08-meta';                   % orbita quasi circolare
%satellite='E08-sp3';
%
% satellite='E18-meta';                 % orbita ellittica
satellite='E18-sp3';
%
% satellite='E30-meta';                   % orbita quasi circolare
%satellite='E30-sp3';
%
% satellite='E14-meta';                 % orbita ellittica
% satellite='E14-sp3';
%
% const=read_par('constant_SI.par');              % carica le costanti SI
% astrody=read_par('astrody_SI.par');             % carica parametri astrofisici    
% par_bw=read_box_wing_par('box_wing_meta.par');  % carica i parametri del boxwing
%
if strcmp(satellite,'E18-meta')
    % GSAT0201 	E18 	Ext01 	27977.6 	0.162 	49.850 	52.521 	56.198 	316.069
    % 2016-11-21 00:00:00 UTC4
    % Reference parameter rates 	Nominal slots 	Extended slots
    % d(RAAN)/dt 	-0.02764398 deg/day 	-0.03986760 deg/day
    % d(Arg. peri)/dt 	0.00000000 deg/day 	0.03383184 deg/day
    % d(Mean Anomaly)/dt 	613.72253566 deg/day 	667.86467481 deg/day
    %
    % True anomaly (?) and mean anomaly (M) are related through the eccentric anomaly and the Kepler’s equation. 
    % The true anomaly can be also solved from the mean anomaly by using a series expansion approach of the so-called 
    % equation of the center. For reference orbit computation the series solution can be truncated in the following terms:
    % Note that for circular orbits (eccentricity =0) both mean and true anomalies are identical.
    % keplerian parameters of the two satellites
    %
    a(1)=27977.6e3;   % semimajoraxes m
    e(1)=0.162;       % eccentricy
    I(1)=49.850/180*pi; % Inclination rad
    OMG(1)=52.521/180*pi;   % longitude ascending node (rad)  
    om(1)=56.198/180*pi;    % pericenter argument (rad)   
    MA(1)=316.069/180*pi;   % mean anomaly (rad)
    dOMG=-0.039867092/180*pi;     %% longitude ascending node variation (rad/day)
    dom=0.034373586/180*pi;
    dMA=667.909221051/180*pi;    % mean anomaly rate

    t_launch=mjuliandate([2014 08 22 12 27 00]);      % launch date of FOC 201 e 202
    m_sat=660.977;
elseif strcmp(satellite,'E18-sp3')
    % GSAT0201 	E18 	Ext01 	27977.6 	0.162 	49.850 	52.521 	56.198 	316.069
    % 2016-11-21 00:00:00 UTC4
    % Reference parameter rates 	Nominal slots 	Extended slots
    % d(RAAN)/dt 	-0.02764398 deg/day 	-0.03986760 deg/day
    % d(Arg. peri)/dt 	0.00000000 deg/day 	0.03383184 deg/day
    % d(Mean Anomaly)/dt 	613.72253566 deg/day 	667.86467481 deg/day
    %
    % True anomaly (?) and mean anomaly (M) are related through the eccentric anomaly and the Kepler’s equation. 
    % The true anomaly can be also solved from the mean anomaly by using a series expansion approach of the so-called 
    % equation of the center. For reference orbit computation the series solution can be truncated in the following terms:
    % Note that for circular orbits (eccentricity =0) both mean and true anomalies are identical.
    % keplerian parameters of the two satellites
    %
    a(1)= 27978099.6571506; % 27977.709e3;   % semimajoraxes m
    e(1)= 0.160352416634393; % 0.160;       % eccentricy
    
    I(1)= 50.369485965388/180*pi; % 50.292/180*pi; % Inclination rad
    OMG(1)= 53.5051781754073/180*pi; % 5.343119967963094e+01/180*pi;   % longitude ascending node (rad)  
    om(1)= 50.1835975672884/180*pi; %5.033083659596059e+01/180*pi;    % pericenter argument (rad)   
    MA(1)= 316.069/180*pi; % 316.069/180*pi; % 329.42/180*pi;   % mean anomaly (rad)
    
    dOMG= -0.04000414086/180*pi; %-4.020339545572749e-02/180*pi;     % longitude ascending node variation (rad/day)
    dom= 0.04910776939/180*pi; % 4.934977316892400e-02/180*pi;       %%% Verificare
    dMA=667.909221051/180*pi; %667.86/180*pi;    % mean anomaly rate
    
    t_launch=mjuliandate([2014 08 22 12 27 00]);      % launch date of FOC 201 e 202
    m_sat=660.977;
elseif strcmp(satellite,'E14-meta')
    % GSAT0202 	E14 	Ext02 	27977.6 	0.162 	49.850 	52.521 	56.198 	136.069
    % 2016-11-21 00:00:00 UTC4
    % Reference parameter rates 	Nominal slots 	Extended slots
    % d(RAAN)/dt 	-0.02764398 deg/day 	-0.03986760 deg/day
    % d(Arg. peri)/dt 	0.00000000 deg/day 	0.03383184 deg/day
    % d(Mean Anomaly)/dt 	613.72253566 deg/day 	667.86467481 deg/day
    %
    % True anomaly (?) and mean anomaly (M) are related through the eccentric anomaly and the Kepler’s equation. 
    % The true anomaly can be also solved from the mean anomaly by using a series expansion approach of the so-called 
    % equation of the center. For reference orbit computation the series solution can be truncated in the following terms:
    % Note that for circular orbits (eccentricity =0) both mean and true anomalies are identical.
    % keplerian parameters of the two satellites
    %
    a(1)=27977.6e3;   % semimajoraxes m
    e(1)=0.162;       % eccentricy
    I(1)=49.850/180*pi; % Inclination rad
    OMG(1)=52.521/180*pi;   % longitude ascending node (rad)  
    om(1)=56.198/180*pi;    % pericenter argument (rad)   
    MA(1)=136.069/180*pi;   % mean anomaly (rad)
    dOMG=-0.039867092/180*pi;     %% longitude ascending node variation (rad/day)
    dom=0.034373586/180*pi;
    dMA=667.909221051/180*pi;    % mean anomaly rate

    t_launch=mjuliandate([2014 08 22 12 27 00]);      % launch date of FOC 201 e 202
    m_sat=660;
elseif strcmp(satellite,'E14-sp3')
    % GSAT0202 	E14 	Ext02 	27977.6 	0.162 	49.850 	52.521 	56.198 	136.069
    % 2016-11-21 00:00:00 UTC4
    % Reference parameter rates 	Nominal slots 	Extended slots
    % d(RAAN)/dt 	-0.02764398 deg/day 	-0.03986760 deg/day
    % d(Arg. peri)/dt 	0.00000000 deg/day 	0.03383184 deg/day
    % d(Mean Anomaly)/dt 	613.72253566 deg/day 	667.86467481 deg/day
    %
    % True anomaly (?) and mean anomaly (M) are related through the eccentric anomaly and the Kepler’s equation. 
    % The true anomaly can be also solved from the mean anomaly by using a series expansion approach of the so-called 
    % equation of the center. For reference orbit computation the series solution can be truncated in the following terms:
    % Note that for circular orbits (eccentricity =0) both mean and true anomalies are identical.
    % keplerian parameters of the two satellites
    %
    a(1)= 27977624.8278597; % 27977.630e3;   % semimajoraxes m
    e(1)= 0.160802344466966; % 0.161;       % eccentricy
    
    I(1)= 50.309056926228/180*pi; % 50.305/180*pi; % Inclination rad
    OMG(1)= 52.4590038528467/180*pi; % 5.246169539556058e+01/180*pi;   % longitude ascending node (rad)  
    om(1)= 52.0863162552391/180*pi; % 5.201252099853627e+01/180*pi;    % pericenter argument (rad)   
    MA(1)= 136.069/180*pi; % 157.24/180*pi;   % mean anomaly (rad)
   
    dOMG= -0.04002923721/180*pi; % -4.004007525344663e-02/180*pi;     %% longitude ascending node variation (rad/day)
    dom= 0.04784293004/180*pi; % 4.871851878194874e-02/180*pi;
    dMA= 667.909221051/180*pi; % 664.4/180*pi;    % mean anomaly rate
    
    t_launch=mjuliandate([2014 08 22 12 27 00]);      % launch date of FOC 201 e 202
    m_sat=660;
elseif strcmp(satellite,'E08-meta')
    % GSAT0208	E08	C07	29599.8	0.0	56.0	197.632	0.0	120.153
    % 0201	660.977	316.89	-13.48	561.92
    % 0202	662.141	311.61	-12.60	562.31
    % d(RAAN)/dt	-0.02764398 deg/day	-0.03986760 deg/day
    % d(Arg. peri)/dt	0.00000000 deg/day	0.03383184 deg/day
    % d(Mean Anomaly)/dt	613.72253566 deg/day	667.86467481 deg/day
    %
    a(1)=29599.8e3;         % semimajoraxes m
    e(1)=0;                 % eccentricy
    I(1)=56/180*pi;         % Inclination rad
    OMG(1)=197.632/180*pi;  % logitude ascending node (rad)
    om(1)=0/180*pi;         % pericenter argument (rad)
%     om(1)=305.1665/180*pi
    MA(1)=120.153/180*pi;   % mean anomaly (rad)
    dOMG=-0.02764398/180*pi;%% logitude ascending node variation (rad/day)
    dom=0/180*pi;
    dMA=613.72253566/180*pi;% mean anomaly rate (rad/day)
    t_launch=mjuliandate([2015 12 17 11 51 00]);      % launch date 
    m_sat=719;
elseif strcmp(satellite,'E08-sp3')
    % GSAT0208	E08	C07	29599.8	0.0	56.0	197.632	0.0	120.153
    % 0201	660.977	316.89	-13.48	561.92
    % 0202	662.141	311.61	-12.60	562.31
    % d(RAAN)/dt	-0.02764398 deg/day	-0.03986760 deg/day
    % d(Arg. peri)/dt	0.00000000 deg/day	0.03383184 deg/day
    % d(Mean Anomaly)/dt	613.72253566 deg/day	667.86467481 deg/day
    a(1)=29600.351e3;         % semimajoraxes m
    e(1)=2.274e-4;           % eccentricy
    I(1)=54.907/180*pi;         % Inclination rad
    OMG(1)=1.972821385041630e+02/180*pi;  % logitude ascending node (rad)
    om(1)= 360.0/180*pi;      % argument of pericenter
    MA(1)=130.77/180*pi;   % mean anomaly (rad)
    dOMG=-2.773706568296141e-02/180*pi;%% logitude ascending node variation (rad/day)
    dom=0.192/180*pi;
    dMA=613.7/180*pi;% mean anomaly rate (rad/day)
    t_launch=mjuliandate([2015 12 17 11 51 00]);      % launch date 
    m_sat=719;
elseif strcmp(satellite,'E30-meta')
    % GSAT0206	E30	A05	29599.8	0.0	56.0	317.632	0.0	0.153
    % 0201	660.977	316.89	-13.48	561.92
    % 0202	662.141	311.61	-12.60	562.31
    % d(RAAN)/dt	-0.02764398 deg/day	-0.03986760 deg/day
    % d(Arg. peri)/dt	0.00000000 deg/day	0.03383184 deg/day
    % d(Mean Anomaly)/dt	613.72253566 deg/day	667.86467481 deg/day
    %
    a(1)=29599.8e3;         % semimajoraxes m
    e(1)=0;                 % eccentricy
    I(1)=56/180*pi;         % Inclination rad
    OMG(1)=317.632/180*pi;  % logitude ascending node (rad)
    om(1)=0/180*pi;         % pericenter argument (rad)
    MA(1)=0.153/180*pi;   % mean anomaly (rad)
    dOMG=-0.02764398/180*pi;%% logitude ascending node variation (rad/day)
    dom=0/180*pi;
    dMA=613.72253466/180*pi;% mean anomaly rate (rad/day)
    t_launch=mjuliandate([2015 09 11 02 08 00]);      % launch date 
    m_sat=719;
elseif strcmp(satellite,'E30-sp3')
    % GSAT0206	E30	A05	29599.8	0.0	56.0	317.632	0.0	0.153
    % 0201	660.977	316.89	-13.48	561.92
    % 0202	662.141	311.61	-12.60	562.31
    % d(RAAN)/dt	-0.02764398 deg/day	-0.03986760 deg/day
    % d(Arg. peri)/dt	0.00000000 deg/day	0.03383184 deg/day
    % d(Mean Anomaly)/dt	613.72253566 deg/day	667.86467481 deg/day
    %
    a(1)=29600.150e3;         % semimajoraxes m
    e(1)=2.384e-4;           % eccentricy
    I(1)=57.088/180*pi;         % Inclination rad
    OMG(1)=3.176868354804990e+02/180*pi;  % logitude ascending node (rad)
    om(1)=3.176868354804990e+02/180*pi;      % eccentricy
    MA(1)=12.96/180*pi;   % mean anomaly (rad)
    dOMG=-2.717538819554530e-02/180*pi;%% logitude ascending node variation (rad/day)
    dom=0.56/180*pi;
    dMA=613.23/180*pi;% mean anomaly rate (rad/day)
    t_launch=mjuliandate([2015 09 11 02 08 00]);      % launch date 
    m_sat=719;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminary operations - calculate time span and parameters
t_ref=mjuliandate([2016 11 21 00 00 00]);         % time when the Keplerian parameters were calculated. (mjd)
t_start=t_launch;                                 % beginning of the simulation period (mjd) %
% t_start=57754;                                  % Bury et al. 2019/2020
%t_start=57388
t_end=t_start+365*8.5; %365*4; %200; %365*4; %183; %365*1;                % end of the simulation period (mjd)
% t_end=t_start+200;                              % Bury et al. 2019/2020
dt=1/24/60/60;                                    % time resolution in the simulation period (day)
dt=1/24/60;
dt=1/24/60/30/2;
dt=1/24/30
dt=1/24/10
dt=1/24/30;
t=t_start:dt:t_end;                           % time point to calculate
t_par=t-t_ref;                                % diffence between the start of the simulation period and the time at which the keplerian parameters wer calculated
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Massa dei satelliti nei diversi intervalli di tempo %%
%
% E18 o GSAT0201. Issue Date: 2014-08-22, Satellite Mass: 660.977 kg
%
if strcmp(satellite,'E18-meta')
    m_sat = 660.977;
elseif strcmp(satellite,'E18-sp3')
    m_sat = 660.977;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% E14 o GSAT0202. Issue Date: 2014-08-22 fino a 2015-07-30 (MJD 57233), Satellite Mass: 662.646 kg
% Issue Date: 2015-07-31, Satellite Mass: 662.141 kg
%
if strcmp(satellite,'E14-meta')
    m_sat = 662.3935;     % Valore medio fra 662.646 e 662.141 kg
elseif strcmp(satellite,'E14-sp3')
    MJD_0 = 57233;
    m_sat(t<=MJD_0)=662.646;
        m_sat(t>MJD_0)= 662.141;  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% E08 o GSAT0208. Issue Date: 2015-12-17 fino a 2022-02-23 (MJD 59633), Satellite Mass: 709.138 kg
% Issue Date: 2022-02-24, Satellite Mass: 709.135 kg
%
if strcmp(satellite,'E08-meta')
    m_sat = 709.1365;     % Valore medio fra 709.138 e 709.135 kg
elseif strcmp(satellite,'E08-sp3')
    MJD_0 = 59633;
    m_sat(t<=MJD_0)=709.138;
        m_sat(t>MJD_0)= 709.135;  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% E30 o GSAT0206. Issue Date: 2015-09-11 fino al 2021-10-28 (MJD 59515), Satellite Mass: 707.740 kg
% Issue Date: 2021-10-29, Satellite Mass: 707.734 kg
%
if strcmp(satellite,'E08-meta')
    m_sat = 707.737;     % Valore medio fra 707.740 e 707.734 kg
elseif strcmp(satellite,'E08-sp3')
    MJD_0 = 59515;
    m_sat(t<=MJD_0)=707.740;
        m_sat(t>MJD_0)= 707.734;  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% plot(t,m_sat.*ones(size(t)))
% figure
% plot(t,m_sat,'o')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the time evolution of keplerian and other parameters at simulation times - ECI frame
OMG_t=OMG+dOMG*t_par;                                                           % values of line node (radiant)
om_t=om+dom*t_par;                                                              % values of perigee (radiant)
MA_t=MA+dMA*t_par;                                                              % values of Mean Anomaly  (radiant)
% MA_t=2*pi+circolare(MA_t);
n_orb=[sin(I)*sin(OMG_t);-sin(I)*cos(OMG_t);cos(I)*ones(1,length(OMG_t))];      % versor normal to the orbit plane in ECI frame
Po=1/sqrt(astrody.GM/a^3)*2*pi                                         % orbital period cirrcola orbit(not used)
Pe=1/dMA*2*pi*86400                                                            % orbital period form Mean Anomaly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Positions of the satellite and the sun in ECI frame
r_sat=satellite_orb_ecc(a,om_t,MA_t,OMG_t,I,e);                                 % satellite position cartesian ECI coordinates vs t
[~,~, r_sun]=sun(t(1,:)+2400000.5);                                             % sun position in cartesian ECI coordinates(Km) vs t
r_sun=r_sun*1000;                                                               % sun position in cartesian ECI coordinates(m) vs t
r=r_sun-r_sat;                                                                  % vector from satellite to sun (m)
dr=sqrt(r(1,:).^2+r(2,:).^2+r(3,:).^2);                                         % Earth Sun distance (m)
nu=umbra(r_sun,r_sat,astrody.RS,astrody.RT);                                    % calculate shadow function 
r_sun=r_sun./(ones(3,1)*sqrt(r_sun(1,:).^2+r_sun(2,:).^2+r_sun(3,:).^2));       % unity vector sun position in cartesian ECI coordinates vs t
r_sat_m=sqrt(r_sat(1,:).^2+r_sat(2,:).^2+r_sat(3,:).^2);                        % distance Earth satellite
r_sat=r_sat./(ones(3,1)*r_sat_m);                                               % unity vector satellite position in cartesian ECI coordinates vs t
r=r./dr;                                                                        % versor from satellite to sun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Positions of the satellite and the sun in orbit frame
r_sun_O=ECI2ORB(OMG_t,om_t,I,r_sun);
r_sat_O=ECI2ORB(OMG_t,om_t,I,r_sat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
u_sun2 =atan2(r_sun_O(2,1:end),r_sun_O(1,1:end));
u_sun =atan2(r_sun_O(2,1:end),r_sun_O(1,1:end));
u_sun(u_sun<0)=u_sun(u_sun<0)+2*pi;

figure
plot(t,u_sun)
title([satellite])
legend('Argument of Latitude of the Sun')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Satellite attitude in ECI frame
% The nominal Galileo spacecrafts attitude is as follows: 
% Z axis towards the Earth Centre (in order to illuminate the Earth with its Navigation Antenna), 
% Y axis is perpendicular to the Sun 
% X axis points towards deep space.  
% In order to maintain the nominal attitude it is necessary to turn (“yaw“) about its Z axis while rotating its solar panels around the Y axis.
% The Yaw angle is the angle between the along track vector X0 and X axis
% Main axes
Z_sat=-r_sat;                  % toward Earth center
Y_sat=cross(r,Z_sat);          % ortogonal to sun position respect to earth
X_sat=cross(Y_sat,Z_sat);      % the third axis is orthogonal to previous two
% 
X0=cross(Z_sat,n_orb);         % along track vector
% and their versors
X_sat=X_sat./sqrt(X_sat(1,:).^2+X_sat(2,:).^2+X_sat(3,:).^2);
Y_sat=Y_sat./sqrt(Y_sat(1,:).^2+Y_sat(2,:).^2+Y_sat(3,:).^2);
Z_sat=Z_sat./sqrt(Z_sat(1,:).^2+Z_sat(2,:).^2+Z_sat(3,:).^2);
X0=X0./sqrt(X0(1,:).^2+X0(2,:).^2+X0(3,:).^2);
% X_sat_old=X_sat;
% Y_sat_old=Y_sat;
%
%
B_sat=cross(r,Y_sat);  % 
%
% Yaw Steering Law
Psi_nc = atan2(dot(r_sun,n_orb),dot(r_sun,cross(r_sat,n_orb)));     % metadata Psi value (in normal condition)
% under some conditons, under the standard law, the yaw angle should change faster
% than possible, in these cases a smoothing to the Yaw steering law must be introduced
% Two auxyliary angles are necessary to define the region of modify Yaw steering law
% beta angle between the orbital plane and the sun, eps angle between
% the satellite position and the projection of the vector Earth Sun on the orbit plane
% beta calculated according patent
% t_eq=mjuliandate([2014 03 20 16 57 49]);       % equinozio di primavera
% beta_pat=asin(sin(23.44/180*pi).*sin(2*pi/365.2564*(t-t_eq)).*cos(I)+cos(2*pi/365.2564*(t-t_eq)).*sin(OMG_t).*sin(I)-cos(23.44/180*pi).*sin(2*pi/365.2564*(t-t_eq)).*cos(OMG_t).*sin(I));
beta=pi/2-angle_2_vect(r_sun,n_orb);         % beta calculated using sun and normal orbit vector
%
x=cross(n_orb,r_sun);                        % first step to calculate eps                   
y=cross(n_orb,x);                            % second step to calculate eps  
eps=acos(dot(r_sat,y));                      % colinearity angle between the sun direction porjection in the orbit plane towrd noon and the satellite
eps2=eps;
eps(eps>pi/2)=pi-eps(eps>pi/2);
%
beta_0=4.1;                                 % 4.1 dai metadati
eps_0=10.0;                                 % 10.0 dai metadati; 10.25 sembra migliore
% Calculate psi behaviour and the  differences respect to standard case.
rr=(abs(beta)<beta_0/180*pi)&(abs(eps)<eps_0/180*pi);                 % condition for modified Psi
rr=find(rr==1);                                                 % indexes of the times of modified conditions
n1=diff(rr);                                                    % identify the starting and ending times of Psi modified    
n2=find(n1~=1);                                                 %   
ind1=[1,n2+1];                                                  % starting times of Psi modified 
ind2=[n2,length(rr)];                                           % ending times of Psi modified 

%
Psi_meta=Psi_nc;                                                % initialize the variable that will contain the Psi value modified according metadatae
delta_Psi_meta=Psi_nc*0;                                        % initialize the variable that will contain the difference between normal and modified Psi
%
if sum(n1>1)>0                                                  % if some periods of Psi must be modified
    ind=[rr(ind1);rr(ind2)];                                    % starting and ending times of Psi modified  
    fc=(abs(eps(ind(1,:)-1))>eps_0/180*pi);                        % last condition in previous epoch
    ind=ind(:,fc);                                              % cancel the period in which eps has reached its condition when beta was in "normal" condition
    for ii=1:size(ind,2)                                        % iteration on the period of modified Psi 
        old=Psi_nc(ind(1,ii):ind(2,ii));                        % original Psi
        Psi_meta(ind(1,ii):ind(2,ii))=90*pi/180*sign(Psi_meta(ind(1,ii)-1))+(Psi_meta(ind(1,ii)-1)-90*pi/180*sign(Psi_meta(ind(1,ii)-1))).*cos(2*pi/5656*[1:ind(2,ii)-ind(1,ii)+1]*dt*86400);       % modified Psi
        delta_Psi_meta(ind(1,ii):ind(2,ii))=old-Psi_meta(ind(1,ii):ind(2,ii));    % difference between original and modified Psi
    end
end


% plot Psi as funcion of time
% figure
% plot(t,Psi_nc/pi*180)              % Psi not corrected
% hold on
% plot(t,Psi_meta/pi*180,'r')        % Psi modified as metadata require
% hold on
% plot(t,beta/pi*180,'k')
% xlabel('Time [MJD]')
% ylabel('\Psi (deg)')
% ylabel('Sun \beta angle [deg]')
% lg={'$\Psi_{nom}$', '$\Psi_{mod}$','$\beta_{\odot}$'} ;
% legend(lg,'Interpreter','latex')
% grid on
% 
% differenza_Psi=Psi_meta-Psi_nc;
% 
% figure
% plot(t,differenza_Psi*180/pi,'m')
% hold on
% plot(t,beta/pi*180,'k')
% xlabel('Time (MJD)')
% ylabel('\Delta\Psi (deg)')
% lg={'$\Psi_{mod} - \Psi_{nom}$','$\beta_{\odot}$'} ;
% legend(lg,'Interpreter','latex')
% 
% figure
% plot(t(1:end-1),diff(Psi_nc)*180/pi)
% xlabel('time (mjd)')
% ylabel('diff Psi_nom (degree)')

%
% Calculate Psi angular velocity
% cutoff=pi;                                      % cutoff=pi*3/4 Massimo
% figure 
% der=abs(diff(Psi_nc));
% der(abs(der)>cutoff)=0;                         % cut off non physical differences due to -pi<Psi<pi
% %der=der/(mean(diff(t)*86400));                  % normalize to delta t
% m_der=max(der)/pi*180/dt;                          % maximum value of not corrected Psi angular velocity 
% plot(t(1:end-1),der/pi*180/dt)                     % plot in degree
% hold on
% der_meta=abs(diff(Psi_meta));
% der_meta(abs(der_meta)>cutoff)=0;               % cut off non physical differences due to -pi<Psi<pi
% % der_meta=der_meta/(mean(diff(t)*86400));        % normalize to delta t
% m_der_meta=max(der_meta/pi*180)/dt;                % maximum value of not corrected Psi metadata angular velocity 
% plot(t(1:end-1),der_meta/pi*180/dt)                % plot in degree
% xlabel('time (mjd)')
% ylabel('Psi (degree/s)')
% lg={['d\Psi/dt (degree/s) (',num2str(m_der),' (degree/s))'],['d\Psi_{meta}/dt (degree/s) (',num2str(m_der_meta),' (degree/s))']} ;
% legend(lg)
% %
% Psi_nc_dot=diff(Psi_nc)*180/pi/dt;                  % Angular speed deg./day
% Psi_meta_dot=diff(Psi_meta)*180/pi/dt;              % Angular speed deg./day
% 
% maxPsincdot=max(Psi_nc_dot)
% minPsincdot=min(Psi_nc_dot)
% 
% maxPsimetadot=max(Psi_meta_dot)
% minPsimetadot=min(Psi_meta_dot)
% 
% media_nc=mean(Psi_nc_dot)
% media_meta=mean(Psi_meta_dot)
% 
% Psi_nc_dot2=diff(Psi_nc_dot)/dt;                    % Angular acceleration deg./d^2
% Psi_meta_dot2=diff(Psi_meta_dot)/dt;                % Angular acceleration deg./d^2
% 
% % Angular speed of Psi
% figure
% plot(t(1:end-1),Psi_nc_dot)
% hold on
% plot(t(1:end-1),Psi_meta_dot,'k')
% xlabel('time (MJD)')
% ylabel('d{\Psi} (degree/d)')
% 
% % Angular acceleration of Psi
% figure
% plot(t(1:end-2),Psi_nc_dot2)
% hold on
% plot(t(1:end-2),Psi_meta_dot2,'k')
% xlabel('time (MJD)')
% ylabel('d2{\Psi} (degree/d^2)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Satellite attitude in ECI frame correction due to modified Psi
% The two rotating axes X_sat and Y_sat are corrected according Psi
% modified expressend in metadata
% X_sat_nc=X_sat;
% Y_sat_nc=Y_sat;
n_pan_o=cross(cross(Y_sat,r),Y_sat);        % versore normale al pannello 
%

%     X1=zeros(3,size(X_sat,2));
%     Y1=zeros(3,size(X_sat,2));
%     for ik=1:size(X_sat,2)
%         mat=[X_sat(1,ik),Y_sat(1,ik),Z_sat(1,ik);X_sat(2,ik),Y_sat(2,ik),Z_sat(2,ik);X_sat(3,ik),Y_sat(3,ik),Z_sat(3,ik)];
%         X1(:,ik)=mat*[cos(delta_Psi_meta(:,ik));-sin(delta_Psi_meta(:,ik));0];
%         Y1(:,ik)=mat*[sin(delta_Psi_meta(:,ik));+cos(delta_Psi_meta(:,ik));0];
%     end
% % 
X_sat_p=X_sat.*cos(delta_Psi_meta)-Y_sat.*sin(delta_Psi_meta);
Y_sat_p=Y_sat.*cos(delta_Psi_meta)+X_sat.*sin(delta_Psi_meta);
% 
X_sat=X_sat_p;
Y_sat=Y_sat_p;
%
% sigma=angle_2_vect(Z_sat,r);            % angle between Z-axis and solar panels normal
% sig=sign(dot(Y_sat,cross(Z_sat,r)));  
% sigma=sigma.*sig;
% ang_pan=pi/2-angle_2_vect(r,Y_sat); % angle between Y-axis and solar panels is should be 90 degree except when Psi is modified
%
n_pan=cross(cross(Y_sat,r),Y_sat);        % versore normale al pannello 
%
% figure
% n=82100:4:82130;
% n=82100:4:82130;
% quiver3(r_sat(1,n),r_sat(2,n),r_sat(3,n),X0(1,n),X0(2,n),X0(3,n),'g')
% hold on
% quiver3(r_sat(1,n),r_sat(2,n),r_sat(3,n),X_sat(1,n),X_sat(2,n),X_sat(3,n),'b')
% % quiver3(r_sat(1,n),r_sat(2,n),r_sat(3,n),Y_sat(1,n),Y_sat(2,n),Y_sat(3,n),'r')
% quiver3(r_sat(1,n),r_sat(2,n),r_sat(3,n),X_sat_old(1,n),X_sat_old(2,n),X_sat_old(3,n),'r')
% quiver3(r_sat(1,n),r_sat(2,n),r_sat(3,n),Y_sat_old(1,n),Y_sat_old(2,n),Y_sat_old(3,n),'r--')
% % quiver3(r_sat(1,n),r_sat(2,n),r_sat(3,n),Z_sat(1,n),Z_sat(2,n),Z_sat(3,n),'m')
% % quiver3(r_sat(1,n),r_sat(2,n),r_sat(3,n),r_sun(1,n),r_sun(2,n),r_sun(3,n),'y')
% quiver3(r_sat(1,n),r_sat(2,n),r_sat(3,n),n_orb(1,n),n_orb(2,n),n_orb(3,n),'black')
% 
% legend('X0','X','Y','Z','sun','n')
% %quiver3(zeros(1,length(10)),zeros(1,length(10)),zeros(1,length(10)),n_orb(1,1:10),n_orb(2,1:10),n_orb(3,1:10),'m')
% 
% figure
% 
% quiver3(r_sat(1,end-10:end),r_sat(2,end-10:end),r_sat(3,end-10:end),X0(1,end-10:end),X0(2,end-10:end),X0(3,end-10:end),'g')
% hold on
% quiver3(r_sat(1,end-10:end),r_sat(2,end-10:end),r_sat(3,end-10:end),X_sat(1,end-10:end),X_sat(2,end-10:end),X_sat(3,end-10:end),'b')
% quiver3(r_sat(1,end-10:end),r_sat(2,end-10:end),r_sat(3,end-10:end),Y_sat(1,end-10:end),Y_sat(2,end-10:end),Y_sat(3,end-10:end),'r')
% quiver3(r_sat(1,end-10:end),r_sat(2,end-10:end),r_sat(3,end-10:end),Z_sat(1,end-10:end),Z_sat(2,end-10:end),Z_sat(3,end-10:end),'m')
% quiver3(r_sat(1,end-10:end),r_sat(2,end-10:end),r_sat(3,end-10:end),r_sun(1,end-10:end),r_sun(2,end-10:end),r_sun(3,end-10:end),'y')
% %quiver3(zeros(1,length(10),1),zeros(1,length(10)),zeros(1,length(10)),n_orb(1,end-10:end),n_orb(2,end-10:end),n_orb(3,end-10:end),'m')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cosine of the angles between the sun and the normal to satellite faces Xp
% Xm Yp Ym Zp Zm SAp Sam
cos_theta_sun.X=dot(r,X_sat);    
cos_theta_sun.Xp=cos_theta_sun.X;
cos_theta_sun.Xm=-cos_theta_sun.X;
cos_theta_sun.Xp(cos_theta_sun.Xp<0)=0;
cos_theta_sun.Xm(cos_theta_sun.Xm<0)=0;
%
cos_theta_sun.Y=dot(r,Y_sat);
cos_theta_sun.Yp=cos_theta_sun.Y;
cos_theta_sun.Ym=-cos_theta_sun.Y;
cos_theta_sun.Yp(cos_theta_sun.Yp<0)=0;
cos_theta_sun.Ym(cos_theta_sun.Ym<0)=0;
%
cos_theta_sun.Z=dot(r,Z_sat);
cos_theta_sun.Zp=cos_theta_sun.Z;
cos_theta_sun.Zm=-cos_theta_sun.Z;
cos_theta_sun.Zp(cos_theta_sun.Zp<0)=0;
cos_theta_sun.Zm(cos_theta_sun.Zm<0)=0;
% %
% % Cosine of the angles between the sun and the axes normal to the panels
cos_theta_sun.SA=dot(r,n_pan); 
cos_theta_sun.SAp=cos_theta_sun.SA;
cos_theta_sun.SAm=-cos_theta_sun.SA;
cos_theta_sun.SAp(cos_theta_sun.SAp<0)=0;
cos_theta_sun.SAm(cos_theta_sun.SAm<0)=0;
%
theta_sun=acos(cos_theta_sun.SAp)*180/pi;


figure
plot(t,cos_theta_sun.SAp,'k')
xlabel('Time (mjd)');
ylabel('$\cos \theta_{\odot}$','Interpreter','latex')
legend('GSAT0201','GSAT0208')
grid on
figure
plot(t,theta_sun,'k')
xlabel('Time (mjd)');
ylabel('$\theta_{\odot}$ (deg)','Interpreter','latex')
legend('GSAT0201','GSAT0208')
grid on



% % Cosine of the angles between the 3 satellite axes and the axes normal to the panels
cos_theta_sun.SA_X=dot(X_sat,n_pan);
cos_theta_sun.SAp_X=cos_theta_sun.SA_X;
cos_theta_sun.SAm_X=-cos_theta_sun.SAp_X;

%
cos_theta_sun.SA_Y=dot(Y_sat,n_pan);
cos_theta_sun.SAp_Y=cos_theta_sun.SA_Y;
cos_theta_sun.SAm_Y=-cos_theta_sun.SAp_Y;

%
cos_theta_sun.SA_Z=dot(Z_sat,n_pan);
cos_theta_sun.SAp_Z=cos_theta_sun.SA_Z;
cos_theta_sun.SAm_Z=-cos_theta_sun.SAp_Z;

%
% s=cos_theta_sun.Xm*1.32+cos_theta_sun.Zm*3.022;
% figure
% plot(t,cos_theta_sun.Xm*1.32)


% Force acting on the different face in body frame
[Fxp,FXxp,FYxp,FZxp]=calc_F_sun('Xp',par_bw,cos_theta_sun,astrody,const,nu,dr);         % Force on X+ face
[Fxm,FXxm,FYxm,FZxm]=calc_F_sun('Xm',par_bw,cos_theta_sun,astrody,const,nu,dr);         % Force on X- face
[Fyp,FXyp,FYyp,FZyp]=calc_F_sun('Yp',par_bw,cos_theta_sun,astrody,const,nu,dr);         % Force on Y+ face
[Fym,FXym,FYym,FZym]=calc_F_sun('Ym',par_bw,cos_theta_sun,astrody,const,nu,dr);         % Force on Y- face
[Fzp,FXzp,FYzp,FZzp]=calc_F_sun('Zp',par_bw,cos_theta_sun,astrody,const,nu,dr);         % Force on Z+ face
[Fzm,FXzm,FYzm,FZzm]=calc_F_sun('Zm',par_bw,cos_theta_sun,astrody,const,nu,dr);         % Force on Z- face
[FXsap,FYsap,FZsap]=calc_F_SA_sun('SAp',par_bw,cos_theta_sun,astrody,const,nu,dr);      % Force on Sa+ face of panel
[FXsam,FYsam,FZsam]=calc_F_SA_sun('SAm',par_bw,cos_theta_sun,astrody,const,nu,dr);      % Force on Sa- face of panel
%
% figure
% plot(FXsap.^2+FYsap.^2+FZsap.^2)
%
Fx_sd=Fxp-Fxm+2*(FXsap+FXsam)+FXxp+FXxm+FXyp+FXym+FXzp+FXzm;  % 
Fy_sd=Fyp-Fym+2*(FYsap+FYsam)+FYxp+FYxm+FYyp+FYym+FYzp+FYzm;
Fz_sd=Fzp-Fzm+2*(FZsap+FZsam)+FZxp+FZxm+FZyp+FZym+FZzp+FZzm;
% 
% Express the forces in gauss reference
Rad_g_sd=-Fz_sd;
At_g_sd=Fx_sd.*cos(Psi_meta)-Fy_sd.*sin(Psi_meta);
Ct_g_sd=-Fx_sd.*sin(Psi_meta)-Fy_sd.*cos(Psi_meta);
%
% % versori gauss
R_g=-Z_sat;
A_g=X_sat.*cos(Psi_meta)-Y_sat.*sin(Psi_meta);
C_g=-X_sat.*sin(Psi_meta)-Y_sat.*cos(Psi_meta);
% 
% figure
% n=82100:4:82130;
% n=82100:4:82130;
% quiver3(r_sat(1,n),r_sat(2,n),r_sat(3,n),R_g(1,n),R_g(2,n),R_g(3,n),'g')
% hold on
% quiver3(r_sat(1,n),r_sat(2,n),r_sat(3,n),A_g(1,n),A_g(2,n),A_g(3,n),'b')
% quiver3(r_sat(1,n),r_sat(2,n),r_sat(3,n),C_g(1,n),C_g(2,n),C_g(3,n),'r')
% quiver3(r_sat(1,n),r_sat(2,n),r_sat(3,n),Z_sat(1,n),Z_sat(2,n),Z_sat(3,n),'m')
% % quiver3(r_sat(1,n),r_sat(2,n),r_sat(3,n),r_sun(1,n),r_sun(2,n),r_sun(3,n),'y')
% %quiver3(r_sat(1,n),r_sat(2,n),r_sat(3,n),n_orb(1,n),n_orb(2,n),n_orb(3,n),'black')
% 
% legend('R_g','A_g','C_g','n')


% Force in sun direction
F_D_sd=(dot(X_sat,r).*Fx_sd+dot(Y_sat,r).*Fy_sd+dot(Z_sat,r).*Fz_sd);  % la direzione D è dal satellite verso il sole e r è dal sole al satellite
F_B_sd=dot(X_sat,B_sat).*Fx_sd+dot(Y_sat,B_sat).*Fy_sd+dot(Z_sat,B_sat).*Fz_sd;
%
% hold on
% plot(t,Psi_meta)

a_D_sd = F_D_sd./m_sat;
a_B_sd = F_B_sd./m_sat;
figure
plot(t,a_B_sd)
legend('acc B')
%%%%%%%


%%%% MEDIE ORA PER ORA PER IL MESE 360x180x30x24
% file2='Alb_GAl_mesi_ora_da_20_anni_E18_multi_5'; % file albedo
% 
% file3='IR_GAl_mesi_ora_da_20_anni_all_R2_new_c'; % file IR
% 
% M1=load(file2);
% Mx_Al=-M1.MM_x; % il vettore originale albedo è calcolato partendo Terra verso il satellite qui più utile satellite terra
% My_Al=-M1.MM_y; 
% Mz_Al=-M1.MM_z; 
% 
% n_r = interp1(M1.r,[1:length(M1.r)],r_sat_m)';  % Given the the orbit ray at each time, calculate in which interval of the caalculated values it follows
% per_r=n_r-floor(n_r);   % calculate the distance from the closer value
% r_1=floor(n_r);           % this is the first value of the interpolation
% r_2=r_1+1;
% 
% %
% M1=load(file3);
% Mx_IR=-M1.MM_x; % il vettore originale IR è calcolato partendo Terra verso il satellite qui più utile satellite terra 
% My_IR=-M1.MM_y; 
% Mz_IR=-M1.MM_z; 
% %
% % Albedo data are given in the not inerzial rotating Earth frame so the
% % time values of the latitude and longitude of the satellite that are in
% % ECI frame are tranformed in the rotating reference.
% %  After the calculation of the albedo values in the rotating frame, the
% %  results must e tranformed again in ECI frame.
% U=mjd2utc(t);
% [LON,LAT,~]=cart2sph(r_sat(1,:),r_sat(2,:),r_sat(3,:)); % posizione angolare  nel riferimento ECI
% LON(LON<0)=LON(LON<0)+2*pi;
% corr_sid=gast(t+2400000.5);   % Dati nel riferimento terrestre rotante
% LON=LON-corr_sid;
% %LON=mod(LON,2*pi);
% LON(LON<0)=LON(LON<0)+2*pi;
% % LON(LON>2*pi)=LON(LON>2*pi)-2*pi;
% 
% 
% d_lon=mean(diff(M1.lon));
% d_lat=mean(diff(M1.lat));
% n_LON=ceil((LON)/d_lon);
% n_LAT=ceil((LAT+pi/2)/d_lat);
% % starting from the matrix lat x lon x month x hour x ray calculated at
% % integer values this parte calcuate the value linearly interpolating the
% % values referred to ray, hour, month. Each value of these variables will
% % be included in an interval of the tabelled value of the matrix. r_1, h_1
% % m_1 are the index in the matrix correspond to the first value to
% % interplate, r_2, h_2 and m_2 to the second value, per_r, per_h, per_m are
% % the relative position in the interval between r_1 and r_2 etc.
% 
% %% fit always on two months, two hours, two rays
% leng_m=eomday(U(:,1),U(:,2));     % lenght of the month
% flag=1-round(U(:,3)./leng_m);         % 0 if the day is in the second part of the month i.e one interpolate with next month
% m1=U(:,2)-flag;    % first month  index
% m1(m1==0)=12;
% m2=m1+1;                            % second month index
% m2(m2==13)=1;   
% length_m=eomday(U(:,1),m1);                     % the data are given at 15 day of the month - the number of days between two different value is the length of the first month 
% per_m=(U(:,3)+length_m.*flag-15)./length_m;     % percentage of month value, due to day if you interpolate or not with previous month
% % interpolation
% per_h=(U(:,5)/60+U(:,6)/3600);              % first hour  index
% h_1=U(:,4)+1;         % hour index
% h_2=h_1+1;          % previous hour index
% h_2(h_2==25)=1; 
% 
% %per_alt=
% 
% 
% % calcolate Albedo interpolating on months, hours and altitude
% % first step - calculate all possibile combination of  past and next 
% % 1 2 -over the month
% % a - b over the hour
% % c d -over the orbit ray
% 
% ind_m1_h1_r1 = sub2ind([180   360    12    24    20] ,n_LAT', n_LON', m1, h_1,r_1);
% ind_m1_h1_r2 = sub2ind([180   360    12    24    20] ,n_LAT', n_LON', m1, h_1,r_2);
% ind_m1_h2_r1 = sub2ind([180   360    12    24    20] ,n_LAT', n_LON', m1, h_2,r_1);
% ind_m1_h2_r2 = sub2ind([180   360    12    24    20] ,n_LAT', n_LON', m1, h_2,r_2);
% ind_m2_h1_r1 = sub2ind([180   360    12    24    20] ,n_LAT', n_LON', m2, h_1,r_1);
% ind_m2_h1_r2 = sub2ind([180   360    12    24    20] ,n_LAT', n_LON', m2, h_1,r_2);
% ind_m2_h2_r1 = sub2ind([180   360    12    24    20] ,n_LAT', n_LON', m2, h_2,r_1);
% ind_m2_h2_r2 = sub2ind([180   360    12    24    20] ,n_LAT', n_LON', m2, h_2,r_2);
% 
% % ind_m1_h1_r1=180*360*12*24*(r_1-1)+180*360*12*(h_1-1)+180*360*(m1-1)+180*(n_LON'-1)+n_LAT';
% % ind_m1_h1_r2=180*360*12*24*(r_2-1)+180*360*12*(h_1-1)+180*360*(m1-1)+180*(n_LON'-1)+n_LAT';
% % ind_m1_h2_r1=180*360*12*24*(r_1-1)+180*360*12*(h_2-1)+180*360*(m1-1)+180*(n_LON'-1)+n_LAT';
% % ind_m1_h2_r2=180*360*12*24*(r_2-1)+180*360*12*(h_2-1)+180*360*(m1-1)+180*(n_LON'-1)+n_LAT';
% % ind_m2_h1_r1=180*360*12*24*(r_1-1)+180*360*12*(h_1-1)+180*360*(m2-1)+180*(n_LON'-1)+n_LAT';
% % ind_m2_h1_r2=180*360*12*24*(r_2-1)+180*360*12*(h_1-1)+180*360*(m2-1)+180*(n_LON'-1)+n_LAT';
% % ind_m2_h2_r1=180*360*12*24*(r_1-1)+180*360*12*(h_2-1)+180*360*(m2-1)+180*(n_LON'-1)+n_LAT';
% % ind_m2_h2_r2=180*360*12*24*(r_2-1)+180*360*12*(h_2-1)+180*360*(m2-1)+180*(n_LON'-1)+n_LAT';
% 
% 
% % calcola indici ore  mese IR
% ind1a=size(Mx_IR,1)*size(Mx_IR,2)*size(Mx_IR,3)*U(:,4)+size(Mx_IR,1)*size(Mx_IR,2)*(U(:,2)-1)+size(Mx_IR,1)*(n_LON'-1)+n_LAT'; % ora 0 è bin 1
% id=U(:,4)+2;
% id(id==25)=1;
% ind2a=size(Mx_IR,1)*size(Mx_IR,2)*size(Mx_IR,3)*(id-1)+size(Mx_IR,1)*size(Mx_IR,2)*(U(:,2)-1)+size(Mx_IR,1)*(n_LON'-1)+n_LAT'; 
% % % calcola indici ore  mese successivo
% UU2=U(:,2);
% UU2=UU2+1;
% UU2(UU2==13)=1;
% % 
% ind1b=size(Mx_IR,1)*size(Mx_IR,2)*size(Mx_IR,3)*U(:,4)+size(Mx_IR,1)*size(Mx_IR,2)*(UU2-1)+size(Mx_IR,1)*(n_LON'-1)+n_LAT'; % ora 0 è bin 1
% ind2b=size(Mx_IR,1)*size(Mx_IR,2)*size(Mx_IR,3)*(id-1)+size(Mx_IR,1)*size(Mx_IR,2)*(UU2-1)+size(Mx_IR,1)*(n_LON'-1)+n_LAT'; 
% %
% 
% 
% % a primo mese - b secondo mese
% % 1 prima ora - 2 ora successiva
% % first step - interpolate ray
% % x component
% Ax_Al_m1_h1=Mx_Al(ind_m1_h1_r1)+(Mx_Al(ind_m1_h1_r2)-Mx_Al(ind_m1_h1_r1)).*per_r;
% Ax_Al_m1_h2=Mx_Al(ind_m1_h2_r1)+(Mx_Al(ind_m1_h2_r2)-Mx_Al(ind_m1_h2_r1)).*per_r;
% Ax_Al_m2_h1=Mx_Al(ind_m2_h1_r1)+(Mx_Al(ind_m2_h1_r2)-Mx_Al(ind_m2_h1_r1)).*per_r;
% Ax_Al_m2_h2=Mx_Al(ind_m2_h2_r1)+(Mx_Al(ind_m2_h2_r2)-Mx_Al(ind_m2_h2_r1)).*per_r;
% % y component
% Ay_Al_m1_h1=My_Al(ind_m1_h1_r1)+(My_Al(ind_m1_h1_r2)-My_Al(ind_m1_h1_r1)).*per_r;
% Ay_Al_m1_h2=My_Al(ind_m1_h2_r1)+(My_Al(ind_m1_h2_r2)-My_Al(ind_m1_h2_r1)).*per_r;
% Ay_Al_m2_h1=My_Al(ind_m2_h1_r1)+(My_Al(ind_m2_h1_r2)-My_Al(ind_m2_h1_r1)).*per_r;
% Ay_Al_m2_h2=My_Al(ind_m2_h2_r1)+(My_Al(ind_m2_h2_r2)-My_Al(ind_m2_h2_r1)).*per_r;
% % z component
% Az_Al_m1_h1=Mz_Al(ind_m1_h1_r1)+(Mz_Al(ind_m1_h1_r2)-Mz_Al(ind_m1_h1_r1)).*per_r;
% Az_Al_m1_h2=Mz_Al(ind_m1_h2_r1)+(Mz_Al(ind_m1_h2_r2)-Mz_Al(ind_m1_h2_r1)).*per_r;
% Az_Al_m2_h1=Mz_Al(ind_m2_h1_r1)+(Mz_Al(ind_m2_h1_r2)-Mz_Al(ind_m2_h1_r1)).*per_r;
% Az_Al_m2_h2=Mz_Al(ind_m2_h2_r1)+(Mz_Al(ind_m2_h2_r2)-Mz_Al(ind_m2_h2_r1)).*per_r;
% % second step - interpolate hour
% % x component
% Ax_Al_m1=Ax_Al_m1_h1+(Ax_Al_m1_h2-Ax_Al_m1_h1).*per_h;
% Ax_Al_m2=Ax_Al_m2_h1+(Ax_Al_m2_h2-Ax_Al_m2_h1).*per_h;
% % y component
% Ay_Al_m1=Ay_Al_m1_h1+(Ay_Al_m1_h2-Ay_Al_m1_h1).*per_h;
% Ay_Al_m2=Ay_Al_m2_h1+(Ay_Al_m2_h2-Ay_Al_m2_h1).*per_h;
% % z component
% Az_Al_m1=Az_Al_m1_h1+(Az_Al_m1_h2-Az_Al_m1_h1).*per_h;
% Az_Al_m2=Az_Al_m2_h1+(Az_Al_m2_h2-Az_Al_m2_h1).*per_h;
% % third step -interpolate month
% Ax_Al=Ax_Al_m1+(Ax_Al_m2-Ax_Al_m1).*per_m;
% Ay_Al=Ay_Al_m1+(Ay_Al_m2-Ay_Al_m1).*per_m;
% Az_Al=Az_Al_m1+(Az_Al_m2-Az_Al_m1).*per_m;
% 
% Ax_IR_1=Mx_IR(ind1a)+(Mx_IR(ind2a)-Mx_IR(ind1a)).*(((U(:,5)/60+U(:,6)/3600)));
% Ay_IR_1=My_IR(ind1a)+(My_IR(ind2a)-My_IR(ind1a)).*(((U(:,5)/60+U(:,6)/3600)));
% Az_IR_1=Mz_IR(ind1a)+(Mz_IR(ind2a)-Mz_IR(ind1a)).*(((U(:,5)/60+U(:,6)/3600)));
% 
% Ax_IR_2=Mx_IR(ind1b)+(Mx_IR(ind2b)-Mx_IR(ind1b)).*(((U(:,5)/60+U(:,6)/3600)));
% Ay_IR_2=My_IR(ind1b)+(My_IR(ind2b)-My_IR(ind1b)).*(((U(:,5)/60+U(:,6)/3600)));
% Az_IR_2=Mz_IR(ind1b)+(Mz_IR(ind2b)-Mz_IR(ind1b)).*(((U(:,5)/60+U(:,6)/3600)));
% 
% Ax_IR=Ax_IR_1+(Ax_IR_2-Ax_IR_1).*(U(:,3)./eomday(U(:,1),U(:,2)));
% Ay_IR=Ay_IR_1+(Ay_IR_2-Ay_IR_1).*(U(:,3)./eomday(U(:,1),U(:,2)));
% Az_IR=Az_IR_1+(Az_IR_2-Az_IR_1).*(U(:,3)./eomday(U(:,1),U(:,2)));
% 
% % Calcolato all'altezza dei Galileo in riferimento terrestre per passare a
% % ECI considero che il riferimento eci è ruotato di -corr_sid rispetto al
% % terrestrre
% % Av=[Ax_Al';Ay_Al';Az_Al'];
% % [LON,LAT,R]=cart2sph(Ax_Al',Ay_Al',Az_Al'); % dati nel riferimento ECI
% % LON=LON+corr_sid;
% % LON(LON>pi)=LON(LON>pi)-2*pi;
% % 
% % [Av(1,:),Av(2,:),Av(3,:)]=sph2cart(LON,LAT,R);
% 
% 
% Av_Al(1,:)=Ax_Al'.*cos(corr_sid)-Ay_Al'.*sin(corr_sid);
% Av_Al(2,:)=Ax_Al'.*sin(corr_sid)+Ay_Al'.*cos(corr_sid);
% Av_Al(3,:)=Az_Al';                
% 
% Av_IR(1,:)=Ax_IR'.*cos(corr_sid)-Ay_IR'.*sin(corr_sid);
% Av_IR(2,:)=Ax_IR'.*sin(corr_sid)+Ay_IR'.*cos(corr_sid);
% Av_IR(3,:)=Az_IR';          
% 
% 
% 
% %
% % projection of the three axes on albedo source direction 
% 
% Av_Al_m=sqrt(Av_Al(1,:).^2+Av_Al(2,:).^2+Av_Al(3,:).^2);
% Av_Al_v=Av_Al./Av_Al_m;
% m=find(Av_Al_m==0);
% Av_Al_v(:,Av_Al_m==0)=zeros(3,length(m));
% 
% 
% Av_IR_m=sqrt(Av_IR(1,:).^2+Av_IR(2,:).^2+Av_IR(3,:).^2);
% Av_IR_v=Av_IR./Av_IR_m;
% m=find(Av_IR_m==0);
% Av_IR_v(:,Av_IR_m==0)=zeros(3,length(m));
% 
% 
% %Av_Al_v=-r_sat;
% Test=dot(r_sat,Av_Al_v);
% 
% % Rad_g_sam=dot(Av_Al,R_g);
% % At_g_sam=dot(Av_Al,A_g);
% % Ct_g_sam=dot(Av_Al,C_g);
% % Rad_g_sam_to=sqrt(Rad_g_sam(1,:).^2+At_g_sam(1,:).^2+Ct_g_sam(1,:).^2);
% 
% cos_theta_Al.X=dot(Av_Al_v,X_sat);
% cos_theta_Al.Xp=cos_theta_Al.X;
% cos_theta_Al.Xm=-cos_theta_Al.X;
% cos_theta_Al.Xp(cos_theta_Al.Xp<0)=0;
% cos_theta_Al.Xm(cos_theta_Al.Xm<0)=0;
% %
% cos_theta_Al.Y=dot(Av_Al_v,Y_sat);
% cos_theta_Al.Yp=cos_theta_Al.Y;
% cos_theta_Al.Ym=-cos_theta_Al.Y;
% cos_theta_Al.Yp(cos_theta_Al.Yp<0)=0;
% cos_theta_Al.Ym(cos_theta_Al.Ym<0)=0;
% %
% cos_theta_Al.Z=dot(Av_Al_v,Z_sat);
% cos_theta_Al.Zp=cos_theta_Al.Z;
% cos_theta_Al.Zm=-cos_theta_Al.Z;
% cos_theta_Al.Zp(cos_theta_Al.Zp<0)=0;
% cos_theta_Al.Zm(cos_theta_Al.Zm<0)=0;
% %
% cos_theta_Al.SA=dot(Av_Al_v,n_pan); % angolo fra vettore da satellite a terra (sorgente albedo) e normale specchio
% cos_theta_Al.SAp=cos_theta_Al.SA;
% cos_theta_Al.SAm=-cos_theta_Al.SA;
% cos_theta_Al.SAp(cos_theta_Al.SAp<0)=0;
% cos_theta_Al.SAm(cos_theta_Al.SAm<0)=0;
% %
% cos_theta_Al.SA_X=dot(X_sat,n_pan);
% cos_theta_Al.SAp_X=cos_theta_Al.SA_X;
% cos_theta_Al.SAm_X=-cos_theta_Al.SAp_X;
% %
% cos_theta_Al.SA_Y=dot(Y_sat,n_pan);
% cos_theta_Al.SAp_Y=cos_theta_Al.SA_Y;
% cos_theta_Al.SAm_Y=-cos_theta_Al.SAp_Y;
% % %
% cos_theta_Al.SA_Z=dot(Z_sat,n_pan);
% cos_theta_Al.SAp_Z=cos_theta_Al.SA_Z;
% cos_theta_Al.SAm_Z=-cos_theta_Al.SAp_Z;


%%%%
% [Fxp,FXxp,FYxp,FZxp]=calc_F_alb('Xp',par_bw,cos_theta_Al,astrody,const,nu,dr,Av_Al_m);             % Force on X+ face
% [Fxm,FXxm,FYxm,FZxm]=calc_F_alb('Xm',par_bw,cos_theta_Al,astrody,const,nu,dr,Av_Al_m);         % Force on X- face
% [Fyp,FXyp,FYyp,FZyp]=calc_F_alb('Yp',par_bw,cos_theta_Al,astrody,const,nu,dr,Av_Al_m);         % Force on Y+ face
% [Fym,FXym,FYym,FZym]=calc_F_alb('Ym',par_bw,cos_theta_Al,astrody,const,nu,dr,Av_Al_m);         % Force on Y- face
% [Fzp,FXzp,FYzp,FZzp]=calc_F_alb('Zp',par_bw,cos_theta_Al,astrody,const,nu,dr,Av_Al_m);         % Force on Z+ face
% [Fzm,FXzm,FYzm,FZzm]=calc_F_alb('Zm',par_bw,cos_theta_Al,astrody,const,nu,dr,Av_Al_m);         % Force on Z- face
% [FXsap,FYsap,FZsap]=calc_F_SA_alb('SAp',par_bw,cos_theta_Al,astrody,const,nu,dr,Av_Al_m);    % Force on Sa+ face of panel
% [FXsam,FYsam,FZsam]=calc_F_SA_alb('SAm',par_bw,cos_theta_Al,astrody,const,nu,dr,Av_Al_m);    % Force on Sa+ face of panel
% %
% Fx_al=Fxp-Fxm+2*(FXsap+FXsam)+FXxp+FXxm+FXyp+FXym+FXzp+FXzm;  % ang_pan difference respect to 90°
% Fy_al=Fyp-Fym+2*(FYsap+FYsam)+FYxp+FYxm+FYyp+FYym+FYzp+FYzm;
% Fz_al=Fzp-Fzm+2*(FZsap+FZsam)+FZxp+FZxm+FZyp+FZym+FZzp+FZzm;
% 
% Express the forces in gauss reference
% Rad_g_al=-Fz_al;
% At_g_al=Fx_al.*cos(Psi_meta)-Fy_al.*sin(Psi_meta);
% Ct_g_al=-Fx_al.*sin(Psi_meta)-Fy_al.*cos(Psi_meta);
% %
% F_D_al=dot(X_sat,r).*Fx_al+dot(Y_sat,r).*Fy_al+dot(Z_sat,r).*Fz_al;  % il verso dal sole al satellite richiede un meno
% F_B_al=dot(X_sat,B_sat).*Fx_al+dot(Y_sat,B_sat).*Fy_al+dot(Z_sat,B_sat).*Fz_al;
% 
% 
% % 
% Mod_al_gauss=sqrt(Rad_g_al.^2+At_g_al.^2+Ct_g_al.^2);
% Mod_al_Eci=sqrt(Fx_al.^2+Fy_al.^2+Fz_al.^2);
%%

% %%% IR
% cos_theta_IR.X=dot(Av_IR_v,X_sat);
% cos_theta_IR.Xp=cos_theta_IR.X;
% cos_theta_IR.Xm=-cos_theta_IR.X;
% cos_theta_IR.Xp(cos_theta_IR.Xp<0)=0;
% cos_theta_IR.Xm(cos_theta_IR.Xm<0)=0;
% %
% cos_theta_IR.Y=dot(Av_IR_v,Y_sat);
% cos_theta_IR.Yp=cos_theta_IR.Y;
% cos_theta_IR.Ym=-cos_theta_IR.Y;
% cos_theta_IR.Yp(cos_theta_IR.Yp<0)=0;
% cos_theta_IR.Ym(cos_theta_IR.Ym<0)=0;
% %
% cos_theta_IR.Z=dot(Av_IR_v,Z_sat);
% cos_theta_IR.Zp=cos_theta_IR.Z;
% cos_theta_IR.Zm=-cos_theta_IR.Z;
% cos_theta_IR.Zp(cos_theta_IR.Zp<0)=0;
% cos_theta_IR.Zm(cos_theta_IR.Zm<0)=0;
% %
% cos_theta_IR.SA=dot(Av_IR_v,n_pan); % angolo fra vettore da satellite a terra (sorgente albedo) e normale specchio
% cos_theta_IR.SAp=cos_theta_IR.SA;
% cos_theta_IR.SAm=-cos_theta_IR.SA;
% cos_theta_IR.SAp(cos_theta_IR.SAp<0)=0;
% cos_theta_IR.SAm(cos_theta_IR.SAm<0)=0;
% %
% cos_theta_IR.SA_X=dot(X_sat,n_pan);
% cos_theta_IR.SAp_X=cos_theta_IR.SA_X;
% cos_theta_IR.SAm_X=-cos_theta_IR.SAp_X;
% % cos_theta_IR.SAp_X(cos_theta_IR.SAp_X<0)=0;
% % cos_theta_IR.SAm_X(cos_theta_IR.SAm_X<0)=0;
% 
% %
% cos_theta_IR.SA_Y=dot(Y_sat,n_pan);
% cos_theta_IR.SAp_Y=cos_theta_IR.SA_Y;
% cos_theta_IR.SAm_Y=-cos_theta_IR.SAp_Y;
% % cos_theta_IR.SAp_Y(cos_theta_IR.SAp_Y<0)=0;
% % cos_theta_IR.SAm_Y(cos_theta_IR.SAm_Y<0)=0;
% %
% cos_theta_IR.SA_Z=dot(Z_sat,n_pan);
% cos_theta_IR.SAp_Z=cos_theta_IR.SA_Z;
% cos_theta_IR.SAm_Z=-cos_theta_IR.SAp_Z;
% % cos_theta_IR.SAp_Z(cos_theta_IR.SAp_Z<0)=0;
% % cos_theta_IR.SAm_Z(cos_theta_IR.SAm_Z<0)=0;
% 
% [Fxp,FXxp,FYxp,FZxp]=calc_F_alb('Xp',par_bw,cos_theta_IR,astrody,const,nu,dr,Av_IR_m);             % Force on X+ face
% [Fxm,FXxm,FYxm,FZxm]=calc_F_alb('Xm',par_bw,cos_theta_IR,astrody,const,nu,dr,Av_IR_m);         % Force on X- face
% [Fyp,FXyp,FYyp,FZyp]=calc_F_alb('Yp',par_bw,cos_theta_IR,astrody,const,nu,dr,Av_IR_m);         % Force on Y+ face
% [Fym,FXym,FYym,FZym]=calc_F_alb('Ym',par_bw,cos_theta_IR,astrody,const,nu,dr,Av_IR_m);         % Force on Y- face
% [Fzp,FXzp,FYzp,FZzp]=calc_F_alb('Zp',par_bw,cos_theta_IR,astrody,const,nu,dr,Av_IR_m);         % Force on Z+ face
% [Fzm,FXzm,FYzm,FZzm]=calc_F_alb('Zm',par_bw,cos_theta_IR,astrody,const,nu,dr,Av_IR_m);         % Force on Z- face
% [FXsap,FYsap,FZsap]=calc_F_SA_alb('SAp',par_bw,cos_theta_IR,astrody,const,nu,dr,Av_IR_m);    % Force on Sa+ face of panel
% [FXsam,FYsam,FZsam]=calc_F_SA_alb('SAm',par_bw,cos_theta_IR,astrody,const,nu,dr,Av_IR_m);    % Force on Sa+ face of panel
% %
% Fx_IR=Fxp-Fxm+2*(FXsap+FXsam)+FXxp+FXxm+FXyp+FXym+FXzp+FXzm;  % ang_pan difference respect to 90°
% Fy_IR=Fyp-Fym+2*(FYsap+FYsam)+FYxp+FYxm+FYyp+FYym+FYzp+FYzm;
% Fz_IR=Fzp-Fzm+2*(FZsap+FZsam)+FZxp+FZxm+FZyp+FZym+FZzp+FZzm;
% % 
% 
% % Express the forces in gauss reference
% Rad_g_IR=-Fz_IR;
% At_g_IR=Fx_IR.*cos(Psi_meta)-Fy_IR.*sin(Psi_meta);
% Ct_g_IR=-Fx_IR.*sin(Psi_meta)-Fy_IR.*cos(Psi_meta);
% %
% F_D_IR=dot(X_sat,r).*Fx_IR+dot(Y_sat,r).*Fy_IR+dot(Z_sat,r).*Fz_IR;  % il verso dal satellite al sole richiede un meno
% F_B_IR=dot(X_sat,B_sat).*Fx_IR+dot(Y_sat,B_sat).*Fy_IR+dot(Z_sat,B_sat).*Fz_IR;
% % 
% Mod_IR_gauss=sqrt(Rad_g_IR.^2+At_g_IR.^2+Ct_g_IR.^2);
% Mod_IR_Eci=sqrt(Fx_IR.^2+Fy_IR.^2+Fz_IR.^2);
% 
% Mod_sd_gauss=sqrt(Rad_g_sd.^2+At_g_sd.^2+Ct_g_sd.^2);
% Mod_sd_Eci=sqrt(Fx_sd.^2+Fy_sd.^2+Fz_sd.^2);
% %
% figure
% plot(t,Mod_sd_gauss./m_sat)
% hold on
% plot(t,Mod_al_gauss./m_sat,'black')
% plot(t,Mod_IR_gauss./m_sat,'red')
% legend('Mod acc _SRP','Mod acc Al ','Mod acc IR ')
% ylabel('acc (m/s^2)')
% xlabel('time (mjd)')
% title('Modulus')
% 
% figure
% subplot(2,1,1)
% plot(t,F_D_sd./m_sat)
% hold on
% plot(t,F_D_al./m_sat,'black')
% plot(t,F_D_IR./m_sat,'r')
% legend('Acc D SRP','Acc D albedo 12x24','Acc D IR 12x24')
% ylabel('acc (m/s^2)')
% %xlabel('time (mjd)')
% title('Acceleration in sun direction')
% subplot(2,1,2)
% plot(t,F_B_sd./m_sat)
% hold on
% plot(t,F_B_al./m_sat,'black')
% plot(t,F_B_IR./m_sat,'r')
% legend('Acc B SRP','Acc B albedo 12x24','Acc B IR 12x24')
% title('Acceleration in B direction')
% ylabel('acc (m/s^2)')
% xlabel('time (mjd)')
% 
% figure
% subplot(3,1,1)
% plot(t,Fx_sd./m_sat)
% hold on
% plot(t,Fx_al./m_sat,'black')
% plot(t,Fx_IR./m_sat,'r')
% legend('Acc X SRP','Acc X albedo 12x24','Acc X IR 12x24')
% ylabel('acc (m/s^2)')
% %xlabel('time (mjd)')
% title('Acceleration in X_B direction')
% subplot(3,1,2)
% plot(t,Fy_sd./m_sat)
% hold on
% plot(t,Fy_al./m_sat,'black')
% plot(t,Fy_IR./m_sat,'r')
% legend('Acc Y SRP','Acc Y albedo 12x24','Acc Y IR 12x24')
% title('Acceleration in Y_B direction')
% ylabel('acc (m/s^2)')
% xlabel('time (mjd)')
% subplot(3,1,3)
% plot(t,Fz_sd./m_sat)
% hold on
% plot(t,Fz_al./m_sat,'black')
% plot(t,Fz_IR./m_sat,'r')
% legend('Acc Z SRP','Acc Z albedo 12x24','Acc Z IR 12x24')
% title('Acceleration in Z_B direction')
% ylabel('acc (m/s^2)')
% xlabel('time (mjd)')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% t1=t-t(1);
% 
% [ff0,Amp0,Fase0]=fftpsd(t1,Psi_meta);
% 
% figure
% plot(1./ff0,Amp0,'r')
% xlabel('Periods (days)','FontSize', 14)
% ylabel('Amplitude  (deg)','FontSize', 14)
% title('FFT \Psi_{meta}')
% grid on
% 
% 
% 
% a_x=Fx_sd./m_sat;
% [ff,Amp,Fase]=fftpsd(t1,a_x);
% 
% figure
% semilogx(ff/86400,Amp,'b')
% xlabel('Frequency (Hz)','FontSize', 14)
% ylabel('Amplitude  (m/s^2)','FontSize', 14)
% title('FFT A_X')
% grid on
% saveas(gcf, 'FFT_Ax', 'fig')
% 
% figure
% plot(1./ff,Amp,'r')
% xlabel('Periods (days)','FontSize', 14)
% ylabel('Amplitude  (m/s^2)','FontSize', 14)
% title('FFT A_X')
% grid on
% saveas(gcf, 'FFT_Ax', 'fig')
% 
% figure
% semilogx(1./ff,Amp,'r')
% xlabel('Periods (days)','FontSize', 14)
% ylabel('Amplitude  (m/s^2)','FontSize', 14)
% title('FFT A_X')
% grid on
% 
% figure
% semilogx(ff,Amp,'b')
% xlabel('Frequencies (1/days)','FontSize', 14)
% ylabel('Amplitude  (m/s^2)','FontSize', 14)
% title('FFT acceleration along X-axis')
% grid on
% 
% 
% a_z=Fz_sd./m_sat;
% [ff1,Amp1,Fase1]=fftpsd(t1,a_z);
% 
% figure
% semilogx(ff1/86400,Amp1,'b')
% xlabel('Frequency (Hz)','FontSize', 14)
% ylabel('Amplitude  (m/s^2)','FontSize', 14)
% title('FFT A_Z')
% grid on
% 
% figure
% plot(1./ff1,Amp1,'r')
% xlabel('Periods (days)','FontSize', 14)
% ylabel('Amplitude  (m/s^2)','FontSize', 14)
% title('FFT A_Z')
% grid on
% saveas(gcf, 'FFT_Az', 'fig')
% 
% figure
% semilogx(1./ff1,Amp1,'r')
% xlabel('Periods (days)','FontSize', 14)
% ylabel('Amplitude  (m/s^2)','FontSize', 14)
% title('FFT A_Z')
% grid on
% 
% forb=1/Pe;
% % figure
% % subplot(2,1,1)
% % semilogx(ff/86400/forb,Amp,'b')
% % xlabel('Frequency (Hz)','FontSize', 12)
% % ylabel('Amplitude  (m/s^2)','FontSize', 12)
% % xlim([9e-6 1e-2])
% % title('FFT acceleration along X-axis')
% % grid on
% % 
% % subplot(2,1,2)
% % semilogx(ff1/86400/forb,Amp1,'b')
% % xlabel('Frequency (Hz)','FontSize', 12)
% % ylabel('Amplitude  (m/s^2)','FontSize', 12)
% % xlim([9e-6 1e-2])
% % title('FFT acceleration along Z-axis')
% % grid on
% 
% figure
% subplot(2,1,1)
% semilogx(ff/86400,Amp,'b')
% xlabel('Frequency (Hz)','FontSize', 12)
% ylabel('Amplitude  (m/s^2)','FontSize', 12)
% xlim([9e-6 1e-2])
% title('FFT acceleration along X-axis')
% grid on
% 
% subplot(2,1,2)
% semilogx(ff1/86400,Amp1,'b')
% xlabel('Frequency (Hz)','FontSize', 12)
% ylabel('Amplitude  (m/s^2)','FontSize', 12)
% xlim([9e-6 1e-2])
% title('FFT acceleration along Z-axis')
% grid on
% 
% figure
% subplot(2,1,1)
% semilogx(1./ff,Amp,'b')
% xlabel('Periods (days)','FontSize', 12)
% ylabel('Amplitude  (m/s^2)','FontSize', 12)
% title('FFT acceleration along X-axis')
% grid on
% 
% subplot(2,1,2)
% semilogx(1./ff1,Amp1,'b')
% xlabel('Periods (days)','FontSize', 12)
% ylabel('Amplitude  (m/s^2)','FontSize', 12)
% title('FFT acceleration along Z-axis')
% grid on
% 
% figure
% subplot(2,1,1)
% semilogx(1./ff/(Pe/86400),Amp,'b')
% xlabel('Number orbital periods','FontSize', 12)
% ylabel('Amplitude  (m/s^2)','FontSize', 12)
% title('FFT acceleration along X-axis')
% grid on
% 
% subplot(2,1,2)
% semilogx(1./ff1/(Pe/86400),Amp1,'b')
% xlabel('Number orbital periods','FontSize', 12)
% ylabel('Amplitude  (m/s^2)','FontSize', 12)
% ylim([0 8e-8])
% title('FFT acceleration along Z-axis')
% grid on
% 
% figure
% semilogx(ff1,Amp1,'b')
% xlabel('Frequencies (1/days)','FontSize', 14)
% ylabel('Amplitude  (m/s^2)','FontSize', 14)
% title('FFT acceleration along Z-axis')
% grid on
% 
% figure
% subplot(2,1,1)
% semilogx(ff,Amp,'b')
% xlabel('Frequencies (1/days)','FontSize', 12)
% ylabel('Amplitude  (m/s^2)','FontSize', 12)
% title('FFT acceleration along X-axis')
% grid on
% 
% subplot(2,1,2)
% semilogx(ff1,Amp1,'b')
% xlabel('Frequencies (1/days)','FontSize', 12)
% ylabel('Amplitude  (m/s^2)','FontSize', 12)
% title('FFT acceleration along Z-axis')
% grid on



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A=r_sun_O;
% B=r_sat_O;
% A(3,:)=zeros(1,length(A));
% B(3,:)=zeros(1,length(A));
% A=A./sqrt(A(1,:).^2+A(2,:).^2);
% B=B./sqrt(B(1,:).^2+B(2,:).^2);
% u=angle_2_vect(A,B);
% cu=dot(A,B);
% su=cross(A,B);
% d=sign(su(3,:));
% 
% su=sqrt(su(1,:).^2+su(2,:).^2+su(3,:).^2);
% u=atan2(su,cu);
% u=u.*d;
% 
% a_D_sd=F_D_sd./m_sat;
% betas=beta/pi*180;
% us=u/pi*180;
% a_B_sd=F_B_sd./m_sat;
% a_D_al=F_D_al./m_sat;
% a_B_al=F_B_al./m_sat;
% a_D_IR=F_D_IR./m_sat;
% a_B_IR=F_B_IR./m_sat;
% 
% figure
% scatter(us,betas,30,Rad_g_sd./m_sat,'filled')
% % caxis([-11.5e-8 -9.5e-8])
% xlabel('\Delta u (degree)')
% ylabel('$\beta_{\odot}$ (degree)','Interpreter','latex')
% colormap jet
% colorbar
% title('SRP radial direction')
% 
% figure
% scatter(us,betas,30,At_g_sd./m_sat,'filled')
% % caxis([-11.5e-8 -9.5e-8])
% xlabel('\Delta u (degree)')
% ylabel('$\beta_{\odot}$ (degree)','Interpreter','latex')
% colormap jet
% colorbar
% title('SRP transverse direction')
% 
% figure
% scatter(us,betas,30,Ct_g_sd./m_sat,'filled')
% % caxis([-11.5e-8 -9.5e-8])
% xlabel('\Delta u (degree)')
% ylabel('$\beta_{\odot}$ (degree)','Interpreter','latex')
% colormap jet
% colorbar
% title('SRP out-of-plane direction')
% 
% figure
% scatter(us,betas,30,a_D_sd,'filled')
% % caxis([-11.5e-8 -9.5e-8])
%  caxis([-13.0e-8 -10.0e-8])
% xlabel('\Delta u (degree)')
% ylabel('$\beta_{\odot}$ (degree)','Interpreter','latex')
% colormap jet
% colorbar
% title('SRP dir D')
% 
% figure
% scatter(us,betas,30,a_B_sd,'filled')
% caxis([-5e-9 5e-9])
% colormap jet
% colorbar
% xlabel('\Delta u (degree)')
% ylabel('$\beta_{\odot}$ (degree)','Interpreter','latex')
% title('SRP dir B')
% 
% figure
% scatter(us,betas,30,a_D_al,'filled')
% caxis([-1e-9 1e-9])
% colormap jet
% colorbar
% xlabel('\Delta u (degree)')
% ylabel('$\beta_{\odot}$ (degree)','Interpreter','latex')
% title('Al dir D (sun)')
% 
% figure
% scatter(us,betas,30,a_B_al,'filled')
% caxis([0 6e-10])
% colormap jet
% colorbar
% xlabel('\Delta u (degree)')
% ylabel('$\beta_{\odot}$ (degree)','Interpreter','latex')
% title('Al dir B ')
% 
% figure
% scatter(us,betas,30,a_D_IR,'filled')
% caxis([-1e-9 1e-9])
% colormap jet
% colorbar
% xlabel('\Delta u (degree)')
% ylabel('$\beta_{\odot}$ (degree)','Interpreter','latex')
% title('IR dir D (sun)')
% 
% figure
% scatter(us,betas,30,a_B_IR,'filled')
% caxis([0 6e-10])
% colormap jet
% colorbar
% xlabel('\Delta u (degree)')
% ylabel('$\beta_{\odot}$ (degree)','Interpreter','latex')
% title('IR dir B ')

% save('fig_bury','u','beta','as')

% figure
% plot(u)
% plot3(beta,u,F_D_sd,'.')

%N=140000;
% N=length(a_D_sd);
% wind='rect';
% 
% nD=find(isnan(a_D_al)==1);
% nB=find(isnan(a_B_al)==1);



% [S_a_D,f]=periodogram((a_D_al),[],N,1/(dt*86400));  % spectrun in day
% [S_a_B,f]=periodogram((a_B_al),[],N,1/(dt*86400));  % spectrun in day 
% pspectrum(x,fs,'spectrogram','TimeResolution',0.0256,'Overlap',86,'Leakage',0.875)
%  %[A_fft,freq]=spectrogram(y,wind,fft_n_over,fft_len,1/Tc);
% f=1/(2*Tc)*linspace(0,1,np/2+1);
% loglog(handles.Op1_freq,1./f,abs(fOp1(1:np/2+1)));
% f = 1/dt*(0:(length(F_sun)/2))/length(F_sun);
%
% fs=1/(dt*86400);
% f = (0:N-1)*(fs/N);     % frequency range
% S_a_D_sd=fft(a_D_sd(1:N))/N;
% S_a_B_sd=fft(a_B_sd(1:N))/N;
% S_a_D_al=fft(a_D_al(1:N))/N;
% S_a_B_al=fft(a_B_al(1:N))/N;
% S_a_D_IR=fft(a_D_IR(1:N))/N;
% S_a_B_IR=fft(a_B_IR(1:N))/N;
% 
% a_At_g_al=At_g_al./m_sat;
% a_Ct_g_al=Ct_g_al./m_sat;
% S_a_At_al=fft(a_At_g_al(1:N))/N;
% S_a_Ct_al=fft(a_Ct_g_al(1:N))/N;
% 
% figure
% subplot(2,1,1)
% plot(1./f/Pe,2*abs(S_a_D_al));
% hold on
% plot(1./f/Pe,2*abs(S_a_D_IR));
% xlim([0 1.5])
% xlabel('Period (orb period)');
% ylabel('acc (m/s^2)')
% legend('a D al','a D IR')
% %title('Con eclissi')
% 
% subplot(2,1,2)
% plot(1./f/Pe,2*abs(S_a_B_al));
% hold on
% plot(1./f/Pe,2*abs(S_a_B_IR));
% xlim([0 1.5])
% xlabel('Period (orb period)');
% ylabel('acc (m/s^2)')
% legend('a B al','a B IR')
% 
% figure
% subplot(2,1,1)
% semilogx(f,2*abs(S_a_D_sd),'k');
%  xlim([9e-6 5e-3])
%  ylim([0 6e-9])
% xlabel('Frequency (Hz)');
% ylabel('Amplitude (m/s^2)')
% % legend('a_D')
% title('FFT acceleration along D direction')
% grid on
% 
% subplot(2,1,2)
% semilogx(f,2*abs(S_a_B_sd),'k');
%  xlim([9e-6 5e-3])
%  ylim([0 6e-9])
% xlabel('Frequency (Hz)');
% ylabel('Amplitude (m/s^2)')
% % legend('a_B')
% title('FFT acceleration along B direction')
% grid on
% 
% figure
% subplot(2,1,1)
% plot(1./f/Pe,2*abs(S_a_D_sd));
% % xlim([0.003 1.5])
% % ylim([0 4e-9])
% xlabel('Number orbital periods');
% ylabel('acc (m/s^2)')
% legend('a_D')
% grid on
% 
% subplot(2,1,2)
% plot(1./f/Pe,2*abs(S_a_B_sd));
% % xlim([0.003 1.5])
% % ylim([0 4e-9])
% xlabel('Number orbital periods');
% ylabel('acc (m/s^2)')
% legend('a_B')
% grid on
% % title('senza eclissi')
% 
% figure
% subplot(2,1,1)
% plot(1./f/Pe,2*abs(S_a_At_al));
% xlim([0 1.5])
% xlabel('Period (orb period)');
% ylabel('acc (m/s^2)')
% legend('a_Al_At')
% 
% subplot(2,1,2)
% plot(1./f/Pe,2*abs(S_a_Ct_al));
% xlim([0 1.5])
% xlabel('Period (orb period)');
% ylabel('acc (m/s^2)')
% legend('a_Al_Ct')
% % title('senza eclissi')

% M_b=sqrt(Fx_sd.^2+Fy_sd.^2+Fz_sd.^2);
% M_gaus=sqrt((Rad_g_sd).^2+(At_g_sd).^2+(Ct_g_sd).^2);
% 
% X0_calc=X_sat.*cos(Psi_meta)-Y_sat.*sin(Psi_meta);
% Y0_calc=-X_sat.*sin(Psi_meta)-Y_sat.*cos(Psi_meta);



% figure 
% subplot(3,1,1)
% plot(t,Fx_sd./m_sat)
% title('a_x')
% subplot(3,1,2)
% plot(t,Fy_sd./m_sat)
% title('a_y')
% subplot(3,1,3)
% plot(t,Fz_sd./m_sat)
% title('a_z')


% figure 
% subplot(3,1,1)
% plot(t,Rad_g_sd./m_sat)
% hold on
% plot(t,Rad_g_al./m_sat,'black') 
% plot(t,Rad_g_IR./m_sat,'r') 
% legend('SRP','Albedo 12x24','IR 12x24')
% ylabel('acc (m/s^2)')
% grid on
% title('Radial acceleration')
% subplot(3,1,2)
% plot(t,At_g_sd./m_sat)
% hold on
% plot(t,At_g_al./m_sat,'black')
% plot(t,At_g_IR./m_sat,'r')
% title('Transverse acceleration')
% legend('SRP','Albedo 12x24','IR 12x24')
% ylabel('acc (m/s^2)')
% grid on
% subplot(3,1,3)
% plot(t,Ct_g_sd./m_sat)
% hold on
% plot(t,Ct_g_al./m_sat,'black')
% plot(t,Ct_g_IR./m_sat,'r')
% title('Out-of-plane acceleration')
% legend('SRP','Albedo 12x24','IR 12x24')
% ylabel('acc (m/s^2)')
% xlabel('time (mjd)')
% grid on
% 
% fs=1/dt/86400;
% fs=1/dt
% %fs=0.01;
% 
% % A_ran1 = 1e-8/sqrt(2)*sqrt(fs)*randn(size(t)); 
% % A_ran2 = 1e-9/sqrt(2)*sqrt(fs)*randn(size(t)); 
% 
% A_ran1 = 1e-9/3*randn(size(t)); 
% A_ran2 = 1e-10/2.5*randn(size(t)); 
% 
% 
% figure
% plot(t,A_ran1)
% hold on
% plot(t,A_ran2)
% 
% Rad1=(Rad_g_sd./m_sat)+A_ran1;
% Rad2=(Rad_g_sd./m_sat)+A_ran2;
% 
% % [pn,fn]=psdnorm(Rad,1000,fs);
% % [pb,fb]=psdblack(Rad,1000,fs);
% [pn,fn]=psdnorm(Rad1,length(Rad1),fs);
% [pb,fb]=psdblack(Rad1,length(Rad1),fs);
% 
% % [pn,fn]=psdnorm(A_ran2,length(Rad1),fs);
% % [pb,fb]=psdblack(A_ran2,length(Rad1),fs);
% 
% figure
% loglog(fn,pn.^0.5,'r')
% xlabel('Frequency [Hz]','FontSize', 14)
% ylabel('Amplitude  (m/s^2/sqrt(Hz))','FontSize', 14)
% legend('SRP: Radial acceleration')
% grid on
% 
% figure
% loglog(fb,pb.^0.5,'k')
% xlabel('Frequencies (Hz)','FontSize', 14)
% ylabel('Amplitude  (m/s^2)','FontSize', 14)
% legend('FFT Radial acceleration')
% grid on
% 
% [pn1,fn1]=psdnorm(Rad2,length(Rad2),fs);
% [pb1,fb1]=psdblack(Rad2,length(Rad2),fs);
% 
% 
% figure
% loglog(fn1,pn1.^0.5,'r')
% xlabel('Frequency [Hz]','FontSize', 14)
% ylabel('Amplitude  (m/s^2/sqrt(Hz))','FontSize', 14)
% legend('SRP: Radial acceleration')
% grid on
% 
% figure
% loglog(fb1,pb1.^0.5,'k')
% xlabel('Frequencies (Hz)','FontSize', 14)
% ylabel('Amplitude  (m/s^2)','FontSize', 14)
% legend('FFT Radial acceleration')
% % grid on


% segnaleX1=[t; MA_t; om_t; OMG_t; Rad_g_sd./m_sat; At_g_sd./m_sat; Ct_g_sd./m_sat];
% fprintf(fid1,'\t %23.19e \t %23.19e \t %23.19e \t %23.19e \t %23.19e \t %23.19e \t %23.19e\n',segnaleX1);
% fclose(fid1);
% 
% segnaleX2=[t; MA_t; om_t; OMG_t; Rad_g_al./m_sat; At_g_al./m_sat; Ct_g_al./m_sat];
% fprintf(fid2,'\t %23.19e \t %23.19e \t %23.19e \t %23.19e \t %23.19e \t %23.19e \t %23.19e\n',segnaleX2);
% fclose(fid2);
% 
% segnaleX3=[t; MA_t; om_t; OMG_t; Rad_g_IR./m_sat; At_g_IR./m_sat; Ct_g_IR./m_sat];
% fprintf(fid3,'\t %23.19e \t %23.19e \t %23.19e \t %23.19e \t %23.19e \t %23.19e \t %23.19e\n',segnaleX3);
% fclose(fid3);
% 
% segnaleX4=[t; MA_t; om_t; a_D_sd; a_B_sd];
% fprintf(fid4,'\t %23.19e \t %23.19e \t %23.19e \t %23.19e \t %23.19e\n',segnaleX4);
% fclose(fid4);
% 
% segnaleX5=[t; Fx_sd./m_sat; Fy_sd./m_sat; Fz_sd./m_sat];
% fprintf(fid5,'\t %23.19e \t %23.19e \t %23.19e \t %23.19e\n',segnaleX5);
% fclose(fid5);
% 
% segnaleX6=[t; Fx_al./m_sat; Fy_al./m_sat; Fz_al./m_sat];
% fprintf(fid6,'\t %23.19e \t %23.19e \t %23.19e \t %23.19e\n',segnaleX6);
% fclose(fid6);
% %
% a=5;



function delta=angle_2_vect(v1,v2)
M_v1=v1./(sqrt(v1(1,:).^2+v1(2,:).^2+v1(3,:).^2));
M_v2=v2./(sqrt(v2(1,:).^2+v2(2,:).^2+v2(3,:).^2));
SIN=cross(M_v1,M_v2);
M_SIN=sqrt(SIN(1,:).^2+SIN(2,:).^2+SIN(3,:).^2);
COS=dot(M_v1,M_v2);
delta=atan2(M_SIN,COS); % angle between the sun and Z sat on four quadrant 
end

function [F,FX,FY,FZ]=calc_F_sun(ax,par_bw,cos_theta,astrody,const,nu,d_sun_sat)
%% This routine calculate the force  due to DSR  acting on the face having ax axes as normal
% par_bw - is the structure containing the optical properties of the box-wing model "Area [m2]	Area [m2]	[alpha] [\rho]	[\delta]" for each face
% nu - is shadow function
% cos_theta - contains the cos of the angle between ax normal to the face and each axes, and ax and the sun direction 
F=0;
FX=0;
FY=0;
FZ=0;
Fp=0;
% Force due to radiation pressure according Farinella Milani expression
for i=1:size(par_bw.(ax),1)
    com=par_bw.(ax)(i,1)*(astrody.CS*(astrody.UA./d_sun_sat).^2.*nu).*cos_theta.(ax)/const.c;
%     F=F-2*((par_bw.(ax)(i,4))/3+par_bw.(ax)(i,3).*cos_theta.(ax)).*com;                         % acceleration along normal to panel
          F=F-2*((par_bw.(ax)(i,4))/3+par_bw.(ax)(i,3).*cos_theta.(ax)+par_bw.(ax)(i,2)/3).*com;                         % acceleration along normal to panel con riemissione termica
    Fp=Fp-com.*(par_bw.(ax)(i,2)+par_bw.(ax)(i,4));                                             % acceleration in the direction of the sun    
end
FX=Fp.*cos_theta.X;             % Project the Force in the direction of the sun on the  3 axes 
FY=Fp.*cos_theta.Y;
FZ=Fp.*cos_theta.Z;
end


function [FX,FY,FZ]=calc_F_SA_sun(ax,par_bw,cos_theta,astrody,const,nu,d_sun_sat)
% This routine calculate the force due to DSR acting on the solar panel having ax axes as normal
% par_bw - is the structure containing the optical properties of the box-wing model "Area [m2]	Area [m2]	[alpha] [\rho]	[\delta]" for each face
% nu - is shadow function
% cos_theta - contains the cos of the angle between ax normal to the face and each axes, and ax and the sun direction 
F=0;
Fp=0;
for i=1:size(par_bw.(ax),1)
%    com=par_bw.(ax)(i,1)*(astrody.CS*(astrody.UA./d_sun_sat).^2.*nu).*cos_theta.(ax).*Av_m/const.c;
   com=par_bw.(ax)(i,1)*(astrody.CS*(astrody.UA./d_sun_sat).^2.*nu).*cos_theta.(ax)/const.c;
   %com=par_bw.(ax)(i,1).*cos_theta.(ax).*Av_m/const.c;
   F=F-2*((par_bw.(ax)(i,4))/3+par_bw.(ax)(i,3).*cos_theta.(ax)).*com;                 % acceleration along normal to panel
%    F=F-2*((par_bw.(ax)(i,4))/3+par_bw.(ax)(i,3).*cos_theta.(ax)+par_bw.(ax)(i,2)/3).*com;                         % acceleration along normal to panel con riemissione termica
   Fp=Fp-com.*(par_bw.(ax)(i,2)+par_bw.(ax)(i,4));                                     % acceleration in the direction of the sun
end
FX=Fp.*cos_theta.X+F.*cos_theta.([ax,'_X']);        %.*cos_theta.X+F.*cos_ang_pan.*sin(-sigma);
FY=Fp.*cos_theta.Y+F.*cos_theta.([ax,'_Y']);        %Fp.*cos_theta.Y+F.*sin_ang_pan;
FZ=Fp.*cos_theta.Z+F.*cos_theta.([ax,'_Z']);        %Fp.*cos_theta.Z+F.*cos_ang_pan.*cos(sigma);

end



function [F,FX,FY,FZ]=calc_F_alb(ax,par_bw,cos_theta,astrody,const,nu,d_sun_sat,Av_m)
F=0;
FX=0;
FY=0;
FZ=0;
Fp=0;
for i=1:size(par_bw.(ax),1)
%    com=par_bw.(ax)(i,1)*(astrody.CS*(astrody.UA./d_sun_sat).^2.*nu).*cos_theta.(ax).*Av_m/const.c;
     com=par_bw.(ax)(i,1).*cos_theta.(ax).*Av_m/const.c;
    F=F-2*((par_bw.(ax)(i,4))/3+par_bw.(ax)(i,3).*cos_theta.(ax)).*com;  % acceleration along normal to panel
    Fp=Fp-com.*(par_bw.(ax)(i,2)+par_bw.(ax)(i,4));                                       % acceleration in the direction of the albedo riection
end
FX=Fp.*cos_theta.X;
FY=Fp.*cos_theta.Y;
FZ=Fp.*cos_theta.Z;
end


function [FX,FY,FZ]=calc_F_SA_alb(ax,par_bw,cos_theta,astrody,const,nu,d_sun_sat,Av_m)
F=0;
Fp=0;
for i=1:size(par_bw.(ax),1)
%    com=par_bw.(ax)(i,1)*(astrody.CS*(astrody.UA./d_sun_sat).^2.*nu).*cos_theta.(ax).*Av_m/const.c;
     com=par_bw.(ax)(i,1).*cos_theta.(ax).*Av_m/const.c;
   F=F-2*((par_bw.(ax)(i,4))/3+par_bw.(ax)(i,3).*cos_theta.(ax)).*com;                 % acceleration along normal to panel
   Fp=Fp-com.*(par_bw.(ax)(i,2)+par_bw.(ax)(i,4));                                     % acceleration in the direction of the albedo versor
end
FX=Fp.*cos_theta.X+F.*cos_theta.([ax,'_X']);         %.*cos_theta.X+F.*cos_ang_pan.*sin(-sigma);
FY=Fp.*cos_theta.Y+F.*cos_theta.([ax,'_Y']);       %Fp.*cos_theta.Y+F.*sin_ang_pan;
FZ=Fp.*cos_theta.Z+F.*cos_theta.([ax,'_Z']);    %Fp.*cos_theta.Z+F.*cos_ang_pan.*cos(sigma);
% FX=F.*cos_ang_pan.*sin(sigma);
% FY=F.*sin_ang_pan;
% FZ=F.*cos_ang_pan.*cos(sigma);
end




function [FX,FY,FZ]=old_calc_F_SA_sun_dir(ax,par_bw,cos_theta,sigma,cos_ang_pan,sin_ang_pan,astrody,const,nu,d_sun_sat)
F=0;
Fp=0;
cos_ang_pan(cos_ang_pan<0)=0;
%sin_ang_pan(cos_ang_pan<0)=0;
for i=1:size(par_bw.(ax),1)
    com=par_bw.(ax)(i,1)*(astrody.CS*(astrody.UA./d_sun_sat).^2.*nu).*cos_ang_pan/const.c;
   F=F-2*((par_bw.(ax)(i,4))/3+par_bw.(ax)(i,3).*cos_ang_pan).*com;  % acceleration along normal to panel
   Fp=Fp-com.*(par_bw.(ax)(i,2)+par_bw.(ax)(i,4));                   % acceleration in the direction of the sun
end
FX=Fp.*cos_theta.X+F.*cos_ang_pan.*sin(-sigma);
FY=Fp.*cos_theta.Y+F.*sin_ang_pan;
FZ=Fp.*cos_theta.Z+F.*cos_ang_pan.*cos(sigma);
% FX=F.*cos_ang_pan.*sin(sigma);
% FY=F.*sin_ang_pan;
% FZ=F.*cos_ang_pan.*cos(sigma);
end
