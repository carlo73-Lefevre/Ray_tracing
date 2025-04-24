function nu=umbra(r_sun,r_sat,R_S,R_E)
% 
% input 
% r_sun - sun position in ECI frame  (3xn) vector
% r_sun - satellite position in ECI frame (3xn) vector
% R_S   - sun ray
% R_E   - Earth ray
% 
% output
% nu - shadow function (1xn)
%
% function from Montenbruc Gill pag 81
% 
r=r_sun-r_sat;                                                  % calcuate relative position
a=asin(R_S./sqrt(r(1,:).^2+r(2,:).^2+r(3,:).^2));               % Sun apparent diameter
b=asin(R_E./sqrt(r_sat(1,:).^2+r_sat(2,:).^2+r_sat(3,:).^2));   % Earth apparent diameter
c1=acos(-dot(r_sat(1:3,:),r(1:3,:))./((sqrt(r_sat(1,:).^2+r_sat(2,:).^2+r_sat(3,:).^2)).*sqrt(r(1,:).^2+r(2,:).^2+r(3,:).^2)));
nu=ones(1,size(r_sun,2));
occ=c1<a+b;                     % occultation (tot+partial)
occ_t=(c1<(b-a));               % total occulation 
occ_p=~occ_t&occ;               % partial occultation
nu(occ_t)=0;
x=(c1(occ_p).^2+a(occ_p).^2-b(occ_p).^2)./(2*c1(occ_p));
y=sqrt(a(occ_p).^2-x.^2);
A=a(occ_p).^2.*acos(x./a(occ_p))+b(occ_p).^2.*acos((c1(occ_p)-x)./b(occ_p))-c1(occ_p).*y;
nu(occ_p)=1-(A./(pi.*a(occ_p).^2));
