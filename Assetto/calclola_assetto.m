
addpath('D:\OneDrive - INAF - Istituto Nazionale di Astrofisica\G4S\Matlab_ray\Funzioni\')

[t,X_sat,Y_sat,Z_sat,cos_theta_sun]   =   Assetto_Galileo_pos_sole('E11',55855,56300,1/1440);

figure(1);
clf;
L =640000;
OO = zeros(length(cos_theta_sun.X(1:5:L)),3);
plot_quiv(OO,[cos_theta_sun.X(1:5:L);cos_theta_sun.Y(1:5:L);cos_theta_sun.Z(1:5:L)]',1);
hold on;
X = [1 0 0];Y = [0 1 0]; Z = [0 0 1];
quiver3([0 0 0],[0 0 0],[0 0 0],X,Y,Z,1.5);
text(1.5,0,0,'X')
text(0,1.5,0,'Y')
text(0,0,1.5,'Z')
%%
step = 1;
S1 = [cos_theta_sun.X(1:step:end-1);cos_theta_sun.Y(1:step:end-1);cos_theta_sun.Z(1:step:end-1)];
S2 = [cos_theta_sun.X(2:step:end);cos_theta_sun.Y(2:step:end);cos_theta_sun.Z(2:step:end)];
%%
step = 1;

SUN = [cos_theta_sun.X(1:step:end);cos_theta_sun.Y(1:step:end);cos_theta_sun.Z(1:step:end)];
Xt = [ones(1,length(SUN)).*[1 0 0]'];
Yt = [ones(1,length(SUN)).*[0 1 0]'];
Zt = [ones(1,length(SUN)).*[0 0 1]'];

figure(2)
clf;
ax(1) =nexttile ;      
 plot((t(1:step:end)-t(1))*24*60,rad2deg(acos(dot(SUN,Xt))),'.-');
title('Angolo sole X')
ax(2) = nexttile ;      
plot((t(1:step:end)-t(1))*24*60,rad2deg(acos(dot(SUN,Yt))),'.-r');
title('Angolo sole Y')
ax(3) = nexttile  ;     
plot((t(1:step:end)-t(1))*24*60,rad2deg(acos(dot(SUN,Zt))),'.-c');
title('Angolo sole Z')
linkaxes(ax,'x')
%%
figure(3)
plot(t(1:step:end-1)-t(1),rad2deg(acos(dot(S2,S1))),'.-')
ylabel('Deg')
xlabel('days')
title(['delta angolo Sun = ',string(step*60),' minuti'])

save('Sun_Body.mat','t','SUN');
