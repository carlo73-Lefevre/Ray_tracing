function VO = ECI2ORB(OMG,om,I,VE)
%To transform a vector (nx3) from Earth Centered Inertial frame to Orbital reference frame
n=size(VE,2);
for ii=1:n
    E_EO=[cos(OMG(ii)).*cos(om(ii))-cos(I).*sin(OMG(ii)).*sin(om(ii)),sin(OMG(ii)).*cos(om(ii))+cos(I).*cos(OMG(ii)).*sin(om(ii)),sin(I).*sin(om(ii))
            -cos(OMG(ii)).*sin(om(ii))-cos(I).*sin(OMG(ii)).*cos(om(ii)),cos(I).*cos(OMG(ii)).*cos(om(ii))-sin(OMG(ii)).*sin(om(ii)),sin(I).*cos(om(ii))
            sin(I).*sin(OMG(ii)),-cos(OMG(ii)).*sin(I), cos(I)];
    VO(:,ii)=E_EO*VE(:,ii);
end
end

