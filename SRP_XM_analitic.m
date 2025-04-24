%  clear;
% clc;
path_fold = 'C:\Users\Carlo\OneDrive - INAF - Istituto Nazionale di Astrofisica\G4S\Matlab_ray\Results\';

Visc = dir([path_fold,'Sun_Body_ACC*']); %'FOC*');
[idv1,idv2] = sort(datenum({Visc.date}'));

filenamev =  Visc(idv2(end)).name;

DAT_visco = load(['C:\Users\Carlo\OneDrive - INAF - Istituto Nazionale di Astrofisica\G4S\Matlab_ray\Results\',filenamev]); 
% DAT_visco = load('C:\Users\Carlo\OneDrive - INAF - Istituto Nazionale di Astrofisica\G4S\Matlab_ray\Results\massimo_g2_Sap.mat');



cd('C:\Users\carlo\OneDrive - INAF - Istituto Nazionale di Astrofisica\G4S\Matlab_ray\Results');
FOC = dir('BOX*'); %'FOC*');
[idx1,idx2] = sort(datenum({FOC.date}'));

filename =  FOC(idx2(end)).name;

SRP = load([path_fold,filename]);


idx = find(DAT_visco.time_utc<=SRP.ACC.time(1,end),1,'last');
DT = ceil(idx/ length(SRP.ACC.time));

delta_tot = length(SRP.SRP_body(:,1));

clear ax
figure
clf
ax(1) = subplot(3,2,1);

 plot(SRP.t-SRP.t(1),SRP.SRP_body(:,1),'-');%SRP.ACC.time,
hold on;grid
plot(DAT_visco.time_utc(1:idx)-DAT_visco.time_utc(1),DAT_visco.SRP_body(1:idx,1),'+');%DAT_visco.t(1,1:idx),
xlabel('days')
legend('FEM SRP X','BW ì AccX')
sgtitle(['Confronto Body Acc FEM vs BOX Wing, Starting Time = ',string(datetime(SRP.t(1), 'ConvertFrom', 'modifiedjuliandate'))])

ax(2) = subplot(3,2,3)
plot(SRP.t-SRP.t(1),SRP.SRP_body(:,2));%SRP.ACC.time,
hold on;grid
plot(DAT_visco.time_utc(1:idx)-DAT_visco.time_utc(1),DAT_visco.SRP_body(1:idx,2))
xlabel('days')
legend('FEM SRP Y')

 ax(3) = subplot(3,2,5)
plot(SRP.t-SRP.t(1),SRP.SRP_body(:,3));%SRP.ACC.time,
hold on;grid
plot(DAT_visco.time_utc(1:idx)-DAT_visco.time_utc(1),DAT_visco.SRP_body(1:idx,3),'.');%DAT_visco.t(1,1:idx),
xlabel('days')
legend('FEM SRP Z','BW AccZ')
linkaxes(ax,'x')

ax(4) =subplot(3,2,2)
plot(SRP.t-SRP.t(1),SRP.SRP_body(:,1)-DAT_visco.SRP_body(1:DT:idx,1));%SRP.ACC.time,SRP.t-SRP.t(1),
grid;
hold on;
% plot(SRP.Acc2_SRP(:,1))
legend('difference FEM SRP X vs BW AccX','Acc 2 ref X')

ax(5) =subplot(3,2,4)
plot(SRP.SRP_body(:,2)-DAT_visco.SRP_body(1:DT:idx,2));%SRP.ACC.time,
grid;
hold on;
% plot(SRP.Acc2_SRP(:,2))
legend('difference FEM SRP Y vs BW Visco AccY','Acc 2 ref Y')

ax(6) =subplot(3,2,6)
plot(SRP.SRP_body(:,3) - DAT_visco.SRP_body(1:DT:idx,3));%SRP.ACC.time,
grid;
hold on;
plot(SRP.Acc2_SRP(:,3))
legend('difference FEM SRP Z vs BW Visco AccZ','Acc 2 ref Z')

% linkaxes(ax,'x')
%%
 load([path_fold,'Sat_',filename]);
%%


passo = +1;
% Definisci il percorso per il salvataggio del video
path_dave = ('C:\Users\carlo\OneDrive - INAF - Istituto Nazionale di Astrofisica\G4S\Matlab_ray\Results\');
v = VideoWriter([path_dave,'G4s_BOX_W_vero3.mp4'],'MPEG-4');
v.FrameRate = 30;
v.Quality = 100; 

open(v);  % Apri il file video
fig = figure(1);
% fig.Visible = 'off';
width = 2080; % Larghezza desiderata
height = 1624; % Altezza desiderata
set(gcf, 'Position', [100, 100, width, height]); % Imposta la posizione e le dimensioni
 fig.Visible = 'on';


% for passo = 1:size(Sat_b,2)
 passo
 clf(fig);  % Clear the figure for the new frame  

ax(1) =subplot(1, 2, 1, 'Parent', fig);  

% p1 = patch('Faces',Sat_b(passo).Faces(ray_b{passo,1}(1,2).face_source,:),'Vertices',Sat_b(passo).Vertex); % facce colpite dai riflessi
% p1.FaceColor = [0.7,0.1,0.1];
% axis equal;hold on;

p0 = patch('Faces',Sat_b(passo).Faces(ray_b{passo,1}(1,2).face_source,:),'Vertices',Sat_b(passo).Vertex); % facce colpite dai riflessi
p0.FaceColor = [0.7,0.1,0.1];
axis equal;hold on;
    
% p2 = patch('Faces',Sat_b(passo).Faces(ray2(passo).face_source,:),'Vertices',Sat_b(passo).Vertex); % facce colpite dai riflessi
% p2.FaceColor = 'y';
[azimuth,elevazione,~] = cart2sph(ray_b{passo,1}(1,2).dir(1,1),ray_b{passo,1}(1,2).dir(1,2),ray_b{passo,1}(1,2).dir(1,3));
% view(-90,-90+rad2deg(acos(ray_b{passo}(1).dir(1,3))))
view(-50,30);

% vettori riflessi che illumineranno
% quiver3(Sat_b(passo).Centers_mesh(ray2(passo).face_from,1) ,...
%         Sat_b(passo).Centers_mesh(ray2(passo).face_from,2) ,...
%         Sat_b(passo).Centers_mesh(ray2(passo).face_from,3),...
%         ray2(passo).dir_from(:,1),ray2(passo).dir_from(:,2),ray2(passo).dir_from(:,3),'g')
camzoom(5)



% % % pause(0.01);
ax(2) =subplot(2, 2, 2, 'Parent', fig);  

plot(SRP.SRP_body(1:passo,1),'.');%SRP.ACC.time,SRP.t(1:passo),
hold on;grid
plot(DAT_visco.SRP_body(1:passo,1),'.');%DAT_visco.t(1,1:idx),DAT_visco.time_utc(1:passo),
ylim([5e-8 1e-7])
yyaxis right
% plot(SRP.t(1:passo),SRP.Acc2_SRP(1:passo,1),'k')
ylim([-8e-10 1e-10])
legend('FEM SRP X','BW Visco AccX','SRP X 2 ref')
sgtitle('Confronto Body Acc FEM vs BOX Wing')

ax(3) =subplot(2, 2, 4, 'Parent', fig);  

plot(SRP.SRP_body(1:passo,3),'.');%SRP.ACC.time,SRP.t(1:passo),
hold on;grid
plot(DAT_visco.SRP_body(1:passo,3),'.');%DAT_visco.t(1,1:idx),DAT_visco.time_utc(1:passo),
ylim([-8e-8 1e-7])
yyaxis right
% plot(SRP.t(1:passo),SRP.Acc2_SRP(1:passo,3),'k')
ylim([-5e-10 5e-10])
legend('FEM SRP z','BW Visco Accz','SRP Z 2 ref')
sgtitle('Confronto Body Acc FEM vs BOX Wing')


    %%% Cattura il frame corrente e scrivilo nel video
    frame = getframe(fig);  % Cattura l'immagine della figura corrente
    writeVideo(v, frame);

% end

close(v);

%% no video

passo = 1;
fig = figure(1);
clf;
% fig.Visible = 'off';

view(319,25);

% % % pause(0.01);
ax(2) =subplot(2, 2, 2, 'Parent', fig);  

plot(passo:passo+100,SRP.SRP_body(passo:passo+100,1),'.');%SRP.ACC.time,
hold on;grid
plot(passo:passo+100,DAT_visco.SRP_body(passo:passo+100,1),'.');%DAT_visco.t(1,1:idx),
ylim([5e-8 1e-7])
yyaxis right
% plot(passo:passo+100,SRP.Acc2_SRP(passo:passo+100,1),'k')
ylim([-8e-10 1e-10])
legend('FEM SRP X','BW Visco AccX','SRP X 2 ref')
sgtitle('Confronto Body Acc FEM vs BOX Wing')

ax(3) =subplot(2, 2, 4, 'Parent', fig);  

plot(passo:passo+100,SRP.SRP_body(passo:passo+100,3),'.');%SRP.ACC.time,
hold on;grid
plot(passo:passo+100,DAT_visco.SRP_body(passo:passo+100,3),'.');%DAT_visco.t(1,1:idx),
ylim([-8e-8 1e-7])
yyaxis right
% plot(passo:passo+100,SRP.Acc2_SRP(passo:passo+100,3),'k')
ylim([-5e-10 5e-10])
legend('FEM SRP z','BW Visco Accz','SRP Z 2 ref')
sgtitle('Confronto Body Acc FEM vs BOX Wing')


%%
% passo = 2107;

width = 2080; % Larghezza desiderata
height = 1300; % Altezza desiderata
set(gcf, 'Position', [100, 60, width, height]); % Imposta la posizione e le dimensioni


ax(1) =subplot(1, 2, 1, 'Parent', fig); 

p1= patch('Faces',Sat_b(passo).Faces,'Vertices',Sat_b(passo).Vertex); % facce colpite dai riflessi
p1.FaceColor ='b';
axis equal;hold on;

p0 = patch('Faces',Sat_b(passo).Faces(ray_b{passo,1}(1,2).face_source,:),'Vertices',Sat_b(passo).Vertex); % facce colpite dai riflessi
p0.FaceColor ='y';
axis equal;hold on;


    
% p2 = patch('Faces',Sat_b(passo).Faces(ray2(passo).face_source,:),'Vertices',Sat_b(passo).Vertex); % facce colpite dai riflessi
% p2.FaceColor = 'c';
[azimuth,elevazione,~] = cart2sph(ray_b{passo,1}(1,2).dir(1,1),ray_b{passo,1}(1,2).dir(1,2),ray_b{passo,1}(1,2).dir(1,3));
% view(-90,-90+rad2deg(acos(ray_b{passo}(1).dir(1,3))))


% vettori riflessi che illumineranno
% quiver3(Sat_b(passo).Centers_mesh(ray2(passo).face_from,1) ,...
%         Sat_b(passo).Centers_mesh(ray2(passo).face_from,2) ,...
%         Sat_b(passo).Centers_mesh(ray2(passo).face_from,3),...
%         ray2(passo).dir_from(:,1),ray2(passo).dir_from(:,2),ray2(passo).dir_from(:,3),'g')
camzoom(1)
