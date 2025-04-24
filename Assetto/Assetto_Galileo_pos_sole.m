function [t,X_sat,Y_sat,Z_sat,cos_theta_sun]    =    Assetto_Galileo_pos_sole(satellite,t_start,t_stop,delta_t)
%
%  Il codice calcola l'assetto del satellite Galileo nel sistema ECI e
%  l'angolo della direzione del sole con gli assi del satellite e con le
%  normali ai vari piani.
%  Es. [t,X_sat,Y_sat,Z_sat,cos_theta_sun]   =   Assetto_Galileo_pos_sole('E11',55855,56300,1/1440);
%
%  Jul, 28th 2023
%
% 
% Variabili in ingresso
% satellite    =    nome del satellite SV ID - Es. 'E11'
% t_start    =    inizio del calcolo in mjd
% t_stop    =    fine del calcolo in mjd
% delta_t    =    passo del calcolo in giorni
%
% VARIABILI IN USCITA
% Coordinate degli assi del satellite in ECI frame
% X_sat    =    
% Y_sat    =   
% Z_sat    =   
%
% Angoli fra assi satellite e sole
% cos_theta_sun.X
% cos_theta_sun.Y
% cos_theta_sun.Z
%
% Angoli fra normale a pannello e sole
% cos_theta_sun.SA
% 
% Angoli fra normale a superfici del box e sole
% cos_theta_sun.Xp
% cos_theta_sun.Yp
% cos_theta_sun.Zp
% cos_theta_sun.Xm
% cos_theta_sun.Ym
% cos_theta_sun.Zm
% 
% Angoli fra normale a superfici delle wings e sole
% cos_theta_sun.SAp
% cos_theta_sun.SAm
% cos_theta_sun.SBp
% cos_theta_sun.SBm
%
%  Jul, 27th 2023
%
%  versione:1.0


addpath('D:\OneDrive - INAF - Istituto Nazionale di Astrofisica\G4S\Matlab_ray\Assetto\Assetto_funzioni\');
%
const       =     read_par('constant_SI.par');              % carica le costanti SI
astrody     =     read_par('astrody_SI.par');             % carica parametri astrofisici    
KP         =         readtable('Kp_metadata.txt');
KP.Properties.RowNames    =    string(KP.SVID);
KP_sel   =    KP(char(satellite),:);
%
a   =   KP_sel.SMA_km*1000;               % semiasse maggiore (m)
e   =   KP_sel.e;                         % eccentricità    
I   =   KP_sel.I_deg/180*pi;              % Inclinazione (rad)
OMG   =   KP_sel.RAAN_deg/180*pi;         % Longitudine del nodo ascendente (rad)
om   =   KP_sel.AP_deg/180*pi;            % Argomento del pericentro (rad)
MA   =   KP_sel.MA_deg/180*pi;            % Anomalia media (rad)
dOMG   =   KP_sel.dRAAN_deg_day/180*pi;   % Rate della longitudine del nodo ascendente (rad/day)
dom   =   KP_sel.dAP_deg_day/180*pi;      % Rate dell'argomento del pericentro (rad/day)
dMA   =   KP_sel.dMA_deg_day/180*pi;      % Rate dell'anomalia media (rad/day)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminary operations - calculate time span and parameters
t_ref      =       mjuliandate([2016 11 21 00 00 00]);     % time when the Keplerian parameters were calculated. (mjd)
t          =       t_start:delta_t:t_stop;                     % time point to calculate in mjd
t_par      =       t-t_ref;                                % diffence between the start of the simulation period and the time at which the keplerian parameters wer calculated
%
% Calculate the time evolution of keplerian and other parameters at simulation times - ECI frame
OMG_t     =    OMG+dOMG*t_par;                                                           % values of line node (radiant)
om_t      =    om+dom*t_par;                                                              % values of perigee (radiant)
MA_t      =    MA+dMA*t_par;                                                              % values of Mean Anomaly  (radiant)
n_orb     =    [sin(I)*sin(OMG_t);-sin(I)*cos(OMG_t);cos(I)*ones(1,length(OMG_t))];      % versor normal to the orbit plane in ECI frame
Po     =   1/sqrt(const.G*astrody.MT/a^3)*2*pi;                                         % orbital period cirrcola orbit(not used)
Pe   =   1/dMA*2*pi*86400;                                                            % orbital period form Mean Anomaly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Positions of the satellite and the sun in ECI frame
r_sat   =   satellite_orb_ecc(a,om_t,MA_t,OMG_t,I,e);                                 % satellite position cartesian ECI coordinates vs t
[~,~, r_sun]   =   sun(t(1,:)+2400000.5);                                             % sun position in cartesian ECI coordinates(Km) vs t
r_sun   =   r_sun*1000;                                                               % sun position in cartesian ECI coordinates(m) vs t
r   =   r_sun-r_sat;                                                                  % vector from satellite to sun (m)
dr   =   sqrt(r(1,:).^2+r(2,:).^2+r(3,:).^2);                                         % Earth Sun distance (m)
nu   =   umbra(r_sun,r_sat,astrody.RS,astrody.RT);                                    % calculate shadow function 
nu   =   1;
r_sun   =   r_sun./(ones(3,1)*sqrt(r_sun(1,:).^2+r_sun(2,:).^2+r_sun(3,:).^2));       % unity vector sun position in cartesian ECI coordinates vs t
r_sat_m   =   sqrt(r_sat(1,:).^2+r_sat(2,:).^2+r_sat(3,:).^2);                        % distance Earth satellite
r_sat   =   r_sat./(ones(3,1)*r_sat_m);                                               % unity vector satellite position in cartesian ECI coordinates vs t
r   =   r./dr;                                                                        % versor from satellite to sun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Positions of the satellite and the sun in orbit frame
r_sun_O   =   ECI2ORB(OMG_t,om_t,I,r_sun);
r_sat_O   =   ECI2ORB(OMG_t,om_t,I,r_sat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Satellite attitude in ECI frame
% The nominal Galileo spacecrafts attitude is as follows: 
% Z axis towards the Earth Centre (in order to illuminate the Earth with its Navigation Antenna), 
% Y axis is perpendicular to the Sun 
% X axis points towards deep space.  
% In order to maintain the nominal attitude it is necessary to turn (“yaw“) about its Z axis while rotating its solar panels around the Y axis.
% The Yaw angle is the angle between the along track vector X0 and X axis
% Main axes
Z_sat   =   -r_sat;                  % toward Earth center
Y_sat   =   cross(r,Z_sat);          % ortogonal to sun position respect to earth
X_sat   =   cross(Y_sat,Z_sat);      % the third axis is orthogonal to previous two
% 
X0   =   cross(Z_sat,n_orb);         % along track vector
% and their versors
X_sat   =   X_sat./sqrt(X_sat(1,:).^2+X_sat(2,:).^2+X_sat(3,:).^2);
Y_sat   =   Y_sat./sqrt(Y_sat(1,:).^2+Y_sat(2,:).^2+Y_sat(3,:).^2);
Z_sat   =   Z_sat./sqrt(Z_sat(1,:).^2+Z_sat(2,:).^2+Z_sat(3,:).^2);
X0   =   X0./sqrt(X0(1,:).^2+X0(2,:).^2+X0(3,:).^2);
%
B_sat   =   cross(r,Y_sat);  % 
%
% Yaw Steering Law
beta   =   pi/2-angle_2_vect(r_sun,n_orb);         % beta calculated using sun and normal orbit vector
%
A   =   r_sun_O;
B   =   r_sat_O;
A(3,:)   =   zeros(1,length(A));
B(3,:)   =   zeros(1,length(A));
A   =   A./sqrt(A(1,:).^2+A(2,:).^2);
B   =   B./sqrt(B(1,:).^2+B(2,:).^2);
eta   =   angle_2_vect(A,B);
cu   =   dot(A,B);
su   =   cross(A,B);
d   =   sign(su(3,:));
su   =   sqrt(su(1,:).^2+su(2,:).^2+su(3,:).^2);
eta   =   atan2(su,cu);
eta   =   eta.*d;
betas   =   beta/pi*180;
etas   =   eta/pi*180;
s_0   =   [-sin(eta).*cos(beta);-sin(beta);-cos(eta).*cos(beta)];
Psi_nc    =    atan2(-s_0(2,:)./sqrt(1-s_0(3,:).^2),-s_0(1,:)./sqrt(1-s_0(3,:).^2));     % metadata Psi value (in normal condition)


x   =   cross(n_orb,r_sun);                        % first step to calculate eps                   
y   =   cross(n_orb,x);                            % second step to calculate eps  
eps   =   acos(dot(r_sat,y));                      % colinearity angle between the sun direction porjection in the orbit plane towrd noon and the satellite
eps(eps>pi/2)   =   pi-eps(eps>pi/2);
%
% % Calculate psi behaviour and the  differences respect to standard case.
rr   =   (abs(s_0(1,:))<15/180*pi)&(abs(s_0(2,:))<2/180*pi);         % condition for modified Psi
rr   =   find(rr==1);                                                 % indexes of the times of modified conditions
n1   =   diff(rr);                                                    % identify the starting and ending times of Psi modified    
n2   =   find(n1 ~=   1);                                                 %   
ind1   =   [1,n2+1];                                                  % starting times of Psi modified 
ind2   =   [n2,length(rr)];                                           %
%
Psi_meta   =   Psi_nc;                                                % initialize the variable that will contain the Psi value modified according metadatae
delta_Psi_meta   =   Psi_nc*0;                                        % initialize the variable that will contain the difference between normal and modified Psi
% %
if sum(n1>1)>0                                                  % if some periods of Psi must be modified
    ind   =   [rr(ind1);rr(ind2)];                                    % starting and ending times of Psi modified  
    fc   =   (abs(eps(ind(1,:)-1))>10/180*pi);                        % last condition in previous epoch
    ind   =   ind(:,fc);                                              % cancel the period in which eps has reached its condition when beta was in "normal" condition
    for ii   =   1:size(ind,2)                                        % iteration on the period of modified Psi 
        old   =   Psi_nc(ind(1,ii):ind(2,ii));                        % original Psi
        c   =   cos(pi*abs(s_0(1,ind(1,ii):ind(2,ii)))/sin(15/180*pi));
        Sx   =   s_0(1,ind(1,ii):ind(2,ii));
        Sy   =   0.5*(1+c)*sin(15/180*pi).*sign(s_0(2,ind(1,ii):ind(2,ii)))+0.5*(1-c).*s_0(2,ind(1,ii):ind(2,ii));
        Sz   =   sqrt(s_0(2,ind(1,ii):ind(2,ii))+s_0(1,ind(1,ii):ind(2,ii)));
        Psi_meta(ind(1,ii):ind(2,ii))   =    atan2(-Sy./sqrt(1-Sz.^2),-Sx./sqrt(1-Sz.^2)); 
        %Psi_meta(ind(1,ii):ind(2,ii))   =   90*pi/180*sign(Psi_meta(ind(1,ii)-1))+(Psi_meta(ind(1,ii)-1)-90*pi/180*sign(Psi_meta(ind(1,ii)-1))).*cos(2*pi/5656*[1:ind(2,ii)-ind(1,ii)+1]*dt*86400);       % modified Psi
        delta_Psi_meta(ind(1,ii):ind(2,ii))   =   old-Psi_meta(ind(1,ii):ind(2,ii));    % difference between original and modified Psi
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Satellite attitude in ECI frame correction due to modified Psi
% The two rotating axes X_sat and Y_sat are corrected according Psi
% modified expressend in metadata
X_sat_nc   =   X_sat;
Y_sat_nc   =   Y_sat;
X_sat   =   X_sat_nc.*cos(delta_Psi_meta)-Y_sat_nc.*sin(delta_Psi_meta);
Y_sat   =   Y_sat_nc.*cos(delta_Psi_meta)+X_sat_nc.*sin(delta_Psi_meta);
%
%
n_pan   =   cross(cross(Y_sat,r),Y_sat);        % versore normale al pannello 
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cosine of the angles between the sun and the normal to satellite faces Xp
% Xm Yp Ym Zp Zm SAp Sam
cos_theta_sun.X   =   dot(r,X_sat);    
cos_theta_sun.Xp   =   cos_theta_sun.X;
cos_theta_sun.Xm   =   -cos_theta_sun.X;
cos_theta_sun.Xp(cos_theta_sun.Xp<0)   =   0;
cos_theta_sun.Xm(cos_theta_sun.Xm<0)   =   0;
%
cos_theta_sun.Y   =   dot(r,Y_sat);
cos_theta_sun.Yp   =   cos_theta_sun.Y;
cos_theta_sun.Ym   =   -cos_theta_sun.Y;
cos_theta_sun.Yp(cos_theta_sun.Yp<0)   =   0;
cos_theta_sun.Ym(cos_theta_sun.Ym<0)   =   0;
%
cos_theta_sun.Z   =   dot(r,Z_sat);
cos_theta_sun.Zp   =   cos_theta_sun.Z;
cos_theta_sun.Zm   =   -cos_theta_sun.Z;
cos_theta_sun.Zp(cos_theta_sun.Zp<0)   =   0;
cos_theta_sun.Zm(cos_theta_sun.Zm<0)   =   0;
%
% % Cosine of the angles between the sun and the axes normal to the panels
cos_theta_sun.SA   =   dot(r,n_pan); 
cos_theta_sun.SAp   =   cos_theta_sun.SA;
cos_theta_sun.SAm   =   -cos_theta_sun.SA;
cos_theta_sun.SAp(cos_theta_sun.SAp<0)   =   0;
cos_theta_sun.SAm(cos_theta_sun.SAm<0)   =   0;
%
% % Cosine of the angles between the 3 satellite axes and the axes normal to the panels
cos_theta_sun.SA_X    =   dot(X_sat,n_pan);
cos_theta_sun.SAp_X   =   cos_theta_sun.SA_X;
cos_theta_sun.SAm_X   =   -cos_theta_sun.SAp_X;
%
cos_theta_sun.SA_Y   =   dot(Y_sat,n_pan);
cos_theta_sun.SAp_Y   =   cos_theta_sun.SA_Y;
cos_theta_sun.SAm_Y   =   -cos_theta_sun.SAp_Y;
%
cos_theta_sun.SA_Z    =   dot(Z_sat,n_pan);
cos_theta_sun.SAp_Z   =   cos_theta_sun.SA_Z;
cos_theta_sun.SAm_Z   =   -cos_theta_sun.SAp_Z;
%






function delta   =   angle_2_vect(v1,v2)
M_v1   =   v1./(sqrt(v1(1,:).^2+v1(2,:).^2+v1(3,:).^2));
M_v2   =   v2./(sqrt(v2(1,:).^2+v2(2,:).^2+v2(3,:).^2));
SIN    =   cross(M_v1,M_v2);
M_SIN  =   sqrt(SIN(1,:).^2+SIN(2,:).^2+SIN(3,:).^2);
COS    =   dot(M_v1,M_v2);
delta  =   atan2(M_SIN,COS); % angle between the sun and Z sat on four quadrant 
return


