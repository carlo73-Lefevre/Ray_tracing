clc;close all;
%FOC_3_SAT Foc3_Ansys_materials_5
pat_res = 'C:\Users\Carlo\OneDrive - INAF - Istituto Nazionale di Astrofisica\G4S\Matlab_ray\Results\';
B1 = load([pat_res, 'Box_WING_SAT_area_diff.mat']);L1 ={'Box_WING_SAT_area_diff'};
% B2 = load([pat_res,'FOC_3_SAT.mat']);L1 = [L1,{'FOC_3_SAT'}]; % mesh coerente
 B3 = load([pat_res,'Foc3_Ansys_materials_5.mat']);L1 = [L1,{'Foc3_Ansys_materials_5'}];
%%
figure(1)
clf;
plot(B1.ACC.time,B1.SRP_body(:,1))
hold on;
% plot(B2.SRP_body(1,:))
plot(B3.ACC.time,B3.SRP_body(:,1))
legend(L1)

figure(2)
clf
X = [1 0 0];Y =[0 1 0]; Z=[0 0 1];
passo = +1;

grafica_riflessioni(Sat_b(passo),ray_b(passo,2))
% quiver3([0 0 0],[0 0 0],[0 0 0],X,Y,Z,2);
quiver3(0,0,0,-sat_sun(1,passo),-sat_sun(2,passo),-sat_sun(3,passo),3);
% quiver3(0,0,0,dot(r_sun_sat(:,step),X_sat(:,step)),dot(r_sun_sat(:,step),Y_sat(:,step)),dot(r_sun_sat(:,step),Z_sat(:,step)),10)
% view(-90,-beta_sun_pan(passo)).
view(-60,17)