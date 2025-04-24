function v = satellite_orb_ecc(a,om,MA,OMG,In,ecc)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
eanom=zeros(1,length(MA));
tanom=zeros(1,length(MA));
% for ii=1:length(MA)
[eanom, tanom]=kepler3_1(MA, ecc);
%     [eanom(ii), tanom(ii)]=kepler3(MA(ii), ecc);
% end
%[eanomb, tanomb]=kepler3_1(MA, ecc);
R=a*(1-ecc^2)./(1+ecc*cos(tanom));
%v=R.*[(cos(eanom+om).*cos(OMG)-sin(eanom+om).*cos(In).*sin(OMG));(cos(eanom+om).*sin(OMG)+sin(eanom+om).*cos(In).*cos(OMG)); (sin(eanom+om).*sin(In))];
%v=R.*[(cos(tanom+om).*cos(OMG)-sin(tanom+om).*cos(In).*sin(OMG));(cos(tanom+om).*sin(OMG)+sin(tanom+om).*cos(In).*cos(OMG)); (sin(tanom+om).*sin(In))];
omg=om+tanom;
v=R.*[(cos(omg).*cos(OMG)-sin(omg).*cos(In).*sin(OMG));(cos(omg).*sin(OMG)+sin(omg).*cos(In).*cos(OMG)); (sin(omg).*sin(In))];
end

