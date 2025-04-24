function [eanom, tanom] = kepler3 (manom, ecc)
   
% solve Kepler's equation for elliptic orbits

% Gooding's two iteration method

% input

%  manom = mean anomaly (radians)
%  ecc   = orbital eccentricity (non-dimensional)

% output

%  eanom = eccentric anomaly (radians)
%  tanom = true anomaly (radians)

% Celestial Computing with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pi2 = 2 * pi;

athird = 1 / 3;

ahalf = 1 / 2;

% convert to original algorithm argument

e1 = 1 - ecc;

if (ecc == 0)
   % circular orbit
   
   tanom = manom;
   eanom = manom;
   return;
end

% modulo mean anomaly

emr = manom - pi2 * fix(manom / pi2);
emr(emr < -pi)= emr(emr < -pi) + pi2;
emr(emr > pi)= emr(emr > pi) - pi2;



ee = emr;
ee(ee < 0)= -ee(ee < 0);
eanom=zeros(1,length(ee));
eanom(ee == 0)= ee(ee == 0) + (manom(ee == 0) - emr(ee == 0));
sta = sqrt(1 - ecc * ecc) * sin(eanom(ee == 0));
cta = cos(eanom(ee == 0)) - ecc;
tanom=zeros(1,length(cta));
tanom(ee == 0) = atan3(sta(ee == 0), cta(ee == 0));

n=find(ee~=0);

e = 1 - e1;

% solve cubic equation

w = dcbsol(e, 2 * e1, 3 * ee);

ee = (ee .* ee + (pi - ee) .* w) / pi;

ee(emr < 0) = -ee(emr < 0);


for i = 1:1:2
    fdd = e * sin(ee);
    fddd = e * cos(ee);
    f = emkep(e1, ee) - emr;
    fd = e1 + 2 * e * sin(ahalf*ee).^ 2;
    tst=(e1 + ee .* ee / 6); %>= 0.25)
   f = (ee - fdd) - emr;
   fd = 1 - fddd;
   f(tst>= 0.25)=emkep(e1, ee) - emr;
   fd(tst>= 0.25)=e1 + 2 * e * sin(ahalf * ee) .^ 2;
    dee = f .* fd ./ (ahalf * f .* fdd - fd .* fd);
    w = fd + ahalf * dee .* (fdd + athird * dee .* fddd);
    fd = fd + dee .* (fdd + ahalf * dee .* fddd);
    ee = ee - (f - dee .* (fd - w)) ./ fd;
end

eanom = ee(ee~=0) + (manom(ee~=0) - emr(ee~=0));

% compute true anomaly

sta = sqrt(1 - ecc * ecc) * sin(eanom);
cta = cos(eanom) - ecc;

tanom(ee~=0) = atan3(sta(ee~=0), cta(ee~=0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = dcbsol (a, b, c)

% solve cubic equation function

if (a == 0 & b == 0 | c == 0)
   y = 0;
else
   bsq = b * b;
   d = sqrt(a) * abs(c);
   d = dcubrt(d + sqrt(b * bsq + d.*d)).^2;
   y = 2 * c ./ (d + b + bsq ./ d);
end

%%%%%%%%%%%%%%%%%%%%%%%

function y = dcubrt (x)

% cube root function

athird = 1 / 3;

yy = abs(x);

if (x == 0)
   tmp = 0;
else
   tmp = yy .^ athird;
   tmp = tmp - athird * (tmp - yy./ tmp .^ 2);
   tmp = sign(x) .* tmp;
end

y = tmp;

%%%%%%%%%%%%%%%%%%%%%%%

function y = emkep (e1, ee)

% solves ee - (1 - e1) * sin(ee)
% when (e1, ee) is close to (1, 0)

y = e1 * ee;
ee2 = -ee .* ee;
term = ee .* (1 - e1);

d = 0;

while (1)
  d = d + 2; 
  term = term .* ee2 / (d * (d + 1));
  y0 = y;
  y = y - term;
  if (y0 == y)
     break;
  end
end





