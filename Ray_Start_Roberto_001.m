% Funzione main del ray tracing
%
% Questo codice implementa un algoritmo di Ray Tracing. Inizialmente vengono
% impostati alcuni path che devono essere seguiti dal codice. In seguito,
% viene importato un file relativo ai materiali utilizzati.
% Successivamente, viene creata una classe "SAT" che rappresenta i
% dati statici del satellite. Il file viene poi caricato e montato ruotando
% i pannelli e calcolando le facce complanari.
% Infine, viene implementata la funzione Ray che esegue il Ray Tracing,
% ovvero il calcolo dei raggi di luce che colpiscono il satellite.
% Questo codice può essere utilizzato anche per generare un video delle varie intersezioni.

clear;clc;
% close all;
tic
% path= 'C:\Users\carlo\OneDrive - INAF - Istituto Nazionale di Astrofisica\G4S\Matlab_ray';
path = 'C:\Users\Carlo\OneDrive - INAF - Istituto Nazionale di Astrofisica\G4S\Matlab_ray\';

addpath(path);
addpath([path,'\Funzioni\'])
addpath([path,'\Dati Caomsol\'])
addpath([path,'\Class\'])


% Import material file
% file parametri nome.dat
% file mesh      nome.txt
% file materiali nome.xml
                
% File utilizzabili per test:
% Sat_FEM_Material
% Box_wing_Massimo Box_wing_fine Box_wing_Z [0 3 -1]
% retro_rif retro_rif_low  [-1 -1 -1];
% bullone
% specchio  [0 1 0]
% Cannon Cannon_fine [0 1 -1];
% Lageos_Carlo Lageos_Carlo_Mas Lageos_Carlo_Mas2
% Box_WING_SAT Box_WING_Wings Box_WING_Wings_center
% FOC_3_SAT
% Box_WING_Wings_center
% Box_WING_Wings_center_002
% Bepi_all

% Box_WING_SAT_003
% Box_WING_Wings_center_003

salva_risultati = 'no';

%%% Start Condictions
Sat_inp.Box    = 'bullone';%'Foc3_Ansys_materials_9';%'Foc3_Ansys_materials_8';%'Foc3_Ansys_materials_8';%'Foc3_Ansys_materials_7';%'Foc3_Ansys_materials_6';%Foc3_Ansys_materials_6 'Foc3_Ansys_materials_5';% 'Box_WING_SAT_area_diff';%'Foc3_Ansys_materials_5';%'Box_WING_SAT_area_diff';%'Foc3_Ansys_materials_5'; %'FOC_3_SAT';%;'Box_WING_SAT_area_diff';%'Box_WING_SAT';%'Box_WING_Wings_center_003';%'Foc3_Ansys_materials_4';%'Foc3_Ansys_materials_3';%Foc3_Ansys_materials_2; 'FOC_SAT_mesh004'; %FOC_SAT_mesh001 'FOC_3_SAT';% 'Box_WING_SAT'; % da qui prende la massa e FI
wing = ['no'];

days      = 0.4;  % tempo totale simulazione
dminutes  = 2; %delta tempo in  in minuti


%%%
Sat  = SAT(); % Classe dati satellite statici
Sat1 = SAT(); % Classe dati satellite statici

% Parametri satellite ordinati per material idx
Sat1  = Read_mesh_MAT_txt(Sat_inp.Box,path);

if strcmp(wing,'yes')
    Sat_inp.Wings = 'Box_WING_Wings_center';
end

if  length(fields(Sat_inp)) == 2
    Sat2 = SAT(); % Classe dati satellite statici
    Sat2  = Read_mesh_MAT_txt(Sat_inp.Wings,path);
end

Sat   = SAT();
ray   = c_Rays();

video = 'off';

if strcmp(video,'on')
    path_dave = ('C:\Users\Carlo\OneDrive - INAF - Istituto Nazionale di Astrofisica\G4S\Matlab_ray\Results\');
    v = VideoWriter([path_dave,'G4s_BOX_W_vero.mp4'],'MPEG-4');
    open(v);
end

%%
%  t_launch = utc2mjd([2014 08 22 12 27 00]));   % effemeridi iniziali di E18
addpath('C:\Users\Carlo\OneDrive - INAF - Istituto Nazionale di Astrofisica\G4S\Box_wing\Confronto Box Wing Surf Surf_Visco\Script di Confronto\Massimo_risultati');
path_save = 'C:\Users\Carlo\OneDrive - INAF - Istituto Nazionale di Astrofisica\G4S\Matlab_ray\Results\';

deltaT = dminutes * 1/24/60;
% days   = 0.1;
% Delta_angle = 3 ;% gradi di risoluzione minima

init_Cond = Generazione_coordinate_Satellite_sole('E18',deltaT,days) % (satellite,dt,day_long)

load([path_save,'Sun_Body_dir.mat']);

% !!! verifica la componenet lungo Y_body del sole, dovrebbe essere non 0
sat_sun = sun_body;  % vettore dal Sat to Sun nel ref Body

% angolo_min  = -2200+1000*cos(2*pi*t*0.0058-1.4); % zona sotto la quale aumentare il dettaglio


[azimuth,elevazione,~] = cart2sph(sat_sun(1,:),sat_sun(2,:),sat_sun(3,:));
azimuth_deg    = rad2deg(azimuth);
elevazione_deg = rad2deg(elevazione);

[azimuth_npan,elevazione_npan,~] = cart2sph(n_pan_body(1,:),n_pan_body(2,:),n_pan_body(3,:));
azimuth_pan_deg    = rad2deg(azimuth_npan);
elevazione_pan_deg = rad2deg(elevazione_npan);

beta_sun_pan = -elevazione_pan_deg;  % ho invertito il segno
alfa_sun_pan = azimuth_pan_deg;
%% Calcola la SRP per ogni direzione di sun_sat

pan_ang = elevazione_pan_deg;

% trova le complanari
if strcmp(wing,'no')
        Sat = trova_complanari(Sat1);
        Sat.Normals_mesh = trova_Normali(Sat);
end

if strcmp(wing,'no')
    Sat.materials_idx = Sat.materials_idx+1;
end


%%
clc;
% calocla le normali del Sat 
    % Sat.Normals_mesh = trova_Normali(Sat);
 
kk = 0;
ii = 1;

% tic



while ii <= length(beta_sun_pan) 


% ruota i pannelli solari e trova le complanari
if strcmp(wing,'yes')
  
         [Sat,Ry] = monta_sat(Sat1,Sat2,beta_sun_pan(ii));
         Sat.Normals_mesh = trova_Normali(Sat);
end



    sun_dir       =  -sat_sun(:,ii)'; % va invertito perchè vogliamo la direzione incidente


    kk = kk+1;
     
    ray_input =   repmat(sun_dir,size(Sat.Centers_mesh,1),1);
    ray1_rif  =   raggi_riflessi(Sat.Normals_mesh,ray_input);
    % 

    [Sat,ray,frame]      = Ray_tracing(Sat,sun_dir,ray1_rif,video);

    

    % Sat_b(ii) = Sat;

    [ACC_S1,ACC_N1,Area] = SRP(Sat,ray(1,1),ray(1,2),ray(1,1).dir(1,:),FI_mod(ii),m_sat);

    % disp(strcat('kk =',string(kk),' ii = ',string(ii)))

    ACC.time(kk)        = t(ii);
    ACC.Area(kk)        = Area;
    ACC.beta(kk)        = beta_sun_pan(ii);
    ACC.alf_x(kk)       = alfa_sun_pan(ii);
    ACC.sun(kk,:)       = ACC_S1;
    ACC.sun_mag(kk)     = norm(ACC_S1);

    ACC.N(kk,:)         = ACC_N1;
    ACC.N_mag(kk,:)     = norm(ACC_N1);

    ACC.SRP_mag(kk)     = norm(ACC_S1 + ACC_N1);
    % 
     Sat_b(kk) = Sat;

    % popola ray_b tenendo conto che ray può andara da 1 a 5
    for jj =1:size(ray,2)
        ray_b(kk,jj) = ray(jj);
    end

    idx_ang = 1;

    ii = ii + idx_ang;


    if strcmp(video,'on')
        writeVideo(v,frame);
    end

ntoc(ii) = toc;
 
     updateString = strcat('Total = ',string(length(beta_sun_pan)),' Current value of ii = ', num2str(ii)); %,' deltatime = ',num2str(toc),' tempo medio = ',num2str(sum(ntoc)/(ii-1)));
     disp(updateString);


end

disp(strcat('Total = ',string(toc),' T medio = ', string(toc/ii)));

%% Seconda riflessione

tic
    for jj = size(ray_b,1):-1:1
         ray2(jj) = reflection_test03_002(ray_b(jj,2),Sat_b(jj),1e-7^2);
    end
toc


%% Terza riflessione

tic
    for jj = size(ray2,2):-1:1
         ray3(jj) = reflection_test03_002(ray2(2),Sat_b(jj),1e-7^2);
    end
toc
%% acc seconda riflessione

[ACC_S2,ACC_N2,Acc2n_tot] = SRP_ref(Sat_b(1,1),[],ray2,-sat_sun',m_sat,FI_mod);

%%
%% Interpolazione
% ACC_int.Nx = interp1(ACC.time,ACC.N(:,1),t,"spline");
% ACC_int.Ny = interp1(ACC.time,ACC.N(:,2),t,"spline");
% ACC_int.Nz = interp1(ACC.time,ACC.N(:,3),t,"spline");
% 
% ACC_int.sunx = interp1(ACC.time,ACC.sun(:,1),t,"spline");
% ACC_int.suny = interp1(ACC.time,ACC.sun(:,2),t,"spline");
% ACC_int.sunz = interp1(ACC.time,ACC.sun(:,3),t,"spline");


% SRP_body = [ACC_int.Nx+ACC_int.sunx;ACC_int.Ny+ACC_int.suny;ACC_int.Nz+ACC_int.sunz;];
SRP_body = [ACC.N(:,1)+ACC.sun(:,1)  ACC.N(:,2)+ACC.sun(:,2)  ACC.N(:,3)+ACC.sun(:,3); ];
% SRP_ECI  = [dot(SRP_body,X_sat);dot(SRP_body,Y_sat);dot(SRP_body,Z_sat)];
% SRP_ECI  = [dot(SRP_body',X_sat(:,1:end-1));dot(SRP_body',Y_sat(:,1:end-1));dot(SRP_body',Z_sat(:,1:end-1))];

% Ottieni la data attuale
currentDate = datestr(now, 'yyyymmdd');
% Ottieni l'ora attuale
currentHour = datestr(now, 'HHMMSS');
% Combina la data e l'ora per formare il nome del file
fileName = sprintf([Sat_inp.Box ,'.mat']);
% Ora puoi salvare i tuoi dati utilizzando questo nome del file
pat_res = 'C:\Users\Carlo\OneDrive - INAF - Istituto Nazionale di Astrofisica\G4S\Matlab_ray\Results\';

% save([pat_res,fileName], 'ACC','SRP_body','Sat_b','ray_b');
ora = datestr(now,'mmmm_dd_yyyy_ HH_MM_SS');
% save([pat_res,'Foc3_Ansys_materials_5_SA_2min_7dd_1ref_SAT_ray',ora,'.mat'],'Sat_b','ray_b','-v7.3','-nocompression');
file_name = join([Sat_inp.Box,'_SA_',string(wing),'_',string(dminutes),'min_',string(days),'dd_1ref_ACC',ora,'.mat']);
filename_senza_spazi = strrep(file_name, ' ', '');  

if strcmp(salva_risultati,'yes')
save([pat_res,char(filename_senza_spazi)],'t','ACC','SRP_body','X_sat','Y_sat','Z_sat','sun_body');
save([pat_res,'Sat_',char(filename_senza_spazi)],'t','Sat_b','ray_b');
end

%%
t0 = '2014 08 22 12 27 00';
t_start_utc = datenum(t0,'yyyy mm dd HH MM SS');
T_all =t_start_utc:deltaT:t_start_utc+days;

% 2010-01-01T00:00:00.000

T_Rob =  datetime(datestr(T_all),'Format','yyyy-MM-dd''T''HH:mm:ss.SSS','TimeZone','UTC');


% save([path_save,'ACC_FOC3_solar_pannel_005.mat'],'t','ACC_int','alfa_sun_pan','Sat_b','ray_b');
descrizione_file = 'Accelerazioni FEM ref ECI';
descrizione_file_quat = 'Assetto FEM nel ref ECI';
nomefile = 'FEM_E08_30d_2min_25_09_24';
desc.Satellite = 'GSAT0208' ;

 % save_data_geodine(desc,T_Rob,SRP_ECI,nomefile,descrizione_file)
  % save_data_geodine_quatrnion(desc,T_Rob,ACC_int,nomefile,descrizione_file_quat)
%%
if strcmp(video,'on')
    close(v);
end



% plot(ACC.beta,ACC.SRP_mag,'.-')


% figure
% n=1:1:25;
% quiver3(sun_body(1,n),sun_body(2,n),sun_body(3,n),R_g(1,n),R_g(2,n),R_g(3,n),'g')

figure(11)
clf
plot(ACC.time,ACC.sun(:,1)+ACC.N(:,1),'-.')
grid;
legend('SRP X')


figure(12)
clf
plot(ACC.time,ACC.sun(:,3)+ACC.N(:,3),'-.')
grid;
legend('SRP Z')
% datetick('x','dd-mm-yy HH:MM:SS','keepticks');
%

%  grafica_riflessioni(Sat_b(1),ray_b(1,2))
%  grafica_riflessioni(Sat_b(1),ray_b(1,3))
%  grafica_riflessioni(Sat_b(1),ray_b(1,4))
%  grafica_riflessioni(Sat_b(1),ray_b(1,5))

%%
X = [1 0 0];Y =[0 1 0]; Z=[0 0 1];
passo = 46;

figure(3);
clf;
% grafica_riflessioni(Sat_b(passo),ray_b(passo,2))
grafica_riflessioni(Sat_b(passo),ray_b(passo,2))
legend(strcat('Elevation = ',string(elevazione_pan_deg(passo))))
% quiver3([0 0 0],[0 0 0],[0 0 0],X,Y,Z,2);
% q3 =quiver3(0,0,0,-sat_sun(1,passo),-sat_sun(2,passo),-sat_sun(3,passo),3);
% quiver3(0,0,0,dot(r_sun_sat(:,step),X_sat(:,step)),dot(r_sun_sat(:,step),Y_sat(:,step)),dot(r_sun_sat(:,step),Z_sat(:,step)),10)
% view(-90,-beta_sun_pan(passo)).

view(-70,270+rad2deg(acos(ray_b(passo,1).dir(1,3))))
% view(-90,beta_sun_pan(step))
% -beta_sun_pan(step)
%%

%  figure(2)
% clf;
%  grafica_riflessioni(Sat_b(1),ray_b(1,2))
%  hold on;
% grafica_riflessioni_source(Sat,ray_b(1,2),ray_b(1,3))
% view(-60,30)
%
% figure(3)
% clf;
% grafica_riflessioni_source(Sat_b(1),ray_b(1,4));
% view(-60,30)
% %%
Sat_materiali(Sat_b(passo))

%
% figure(3)
% clf
% grafica_riflessioni(Sat,ray_b(1,1))
%
% figure(4)
% clf;
% grafica_riflessioni(Sat_b(1),ray_b(1,2));
%
%
% figure(5)
% clf
% grafica_riflessioni_source(Sat,ray_b(1,2),ray_b(1,3))
%
% figure(6)
% clf
% grafica_riflessioni_source(Sat,ray_b(1,3),ray_b(1,4))
%
% figure(7)
% clf
% grafica_riflessioni_source(Sat,ray_b(1,4),ray_b(1,5))
%


% figure(3)
% clf;
% p2 = grafica_riflessioni(Sat_b(1),ray_b(1,4));
% p2.FaceColor = 'r';
% view(-60,30)
%
% figure(4)
% clf;
% p2 = grafica_riflessioni(Sat_b(1),ray_b(1,5));
% p2.FaceColor = 'm';
% view(-60,30)


% plot_quiv([0 0 0],[1 0 0],5) % X+ mai al sole
% plot_quiv([0 0 0],[0 0 1],5) % Z+


% save('FOC_centered_003_fine.mat','ACC','Sat_b','ray_b','-mat','-v7.3')


%
% figure(1)
% clf
% surf(ACC.beta,ACC.alf_x,ACC.sun_mag)
% plot(beta_y,ACC.sun_mag,'b--')
% hold on; grid on;
% plot(beta_y,ACC.N_mag,'r--')
% plot(beta_y,ACC.SRP_mag,'g--')
% legend('Sun','Norm','Tot')


toc

%%%% Elenco FUnzioni %%%%%

%% colora facce materiali
% %
function Sat_materiali(Sat)

figure(2)
clf
hold on;
axis equal
box on;

p2 = patch('Faces',Sat.Faces,'Vertices',Sat.Vertex,'EdgeColor','k'); % facce colpite dai riflessi
p2.FaceColor = [0.7,0.7,0.7];
hold on

cmap =[1 0.1 0.1;...
    0.7 0.7 0.7;...
    0 0 1;...
    1 1 0;...
    0 0.5 0.5;...
    0.9 0.8 0.8;...
    1 1 1;...
    0.1 0.5 1;...
    0.2 0.5 0.1;...
    0 0 0;];

L ={};
for ii = 1:max(Sat.materials_idx)
    idx_mat_log = (Sat.materials_idx == ii);  % vertici che contengono il mat
    [idx_row,idx_col] = find(idx_mat_log == 1);

    idx_row_Fn = [];
    for jj=1:length(idx_row)

        [idx_row_F,idx_col_F] = find(Sat.Faces == idx_row(jj));
        idx_row_Fn = [idx_row_Fn;idx_row_F];


    end

    p5 = patch('Faces',Sat.Faces(idx_row_Fn,:),'Vertices',Sat.Vertex,'EdgeColor','none'); % faccia in esame per capire se sta dietro
    p5.FaceColor = cmap(ii,:);
    p5.EdgeLighting =  'gouraud';
    axis equal

    h1(ii) = p5 ;

    L{ii,1} = {(Sat.optical_par(ii).name)};

end

legend(h1,string(L))
title('Satellite Materials','Interpreter','latex')
end

function [Sat,Ry] = monta_sat(Sat1,Sat2,beta_y)
% funzione che unisce i 2 3D, ruota i pannelli solari di un angolo beta_y
% intorno a Y
% e calcola quale faccesono complanari

% ruota pannelli di un angolo alfa
angy = deg2rad(-beta_y+90);


Ry = [cos(angy)  0  sin(angy);...
      0          1          0;...
     -sin(angy)  0  cos(angy)];


% ruota i pannelli solari
% Sat = Sat2;

Sat2.Vertex = [Ry*Sat2.Vertex']';

% Sat1.Vertex = [Ry*Sat1.Vertex']';

% unisci i due body
 Sat = unisci_bodys(Sat1,Sat2);

% determina i centri Mesh (da fare dopo la rotazione=)
Sat               = centri_mesh(Sat);

F_comp            = triangoli_complanari(Sat);
Sat.Face_coplanar = F_comp;

% determina quali facce sono complanari
% L'induce riga è relativo all'indice della Facce in esame che sono
% complanari se il valore è uguale
% idx_faces = [false(length(Sat.Faces),1)]';
% App = (1:length(Sat.Faces))';
% App = cerca_piani_app(idx_faces,Sat.Centri_mesh,App,1,Sat);
% Sat.Face_rige_compl = App;

% Sat = trova_complanari(Sat);

end

function Sat = trova_complanari(Sat)

% determina quali facce sono complanari
% L'induce riga è relativo all'indice della Facce in esame che sono
% complanari se il valore è uguale

%  trova i centri delle mesh
Sat       = centri_mesh(Sat);

Fcomp =  triangoli_complanari(Sat);   % trova le facce complanari e metti un indice comune

Sat.Face_coplanar = Fcomp;

end

function save_data_geodine(desc,t,SRP_ECI,nomefile,Note)
% # Filename:    <sat>_<modello_acc>.acc_est.asc
% # Filetype:    ASCII
% # Author:      <nome_autore_file> (email:<email_autore_file>)
% # Last update: <data_aggiornamento> # Preferibilmente nel formato 'ddd dd mmm yyyy, hh.mm.ss, CET' dato dal comando Unix 'date'
% #
% # Description: Geodyn external accelerations
% #              <descrizione_contenuto_file>
% #              ...
% #
% # Format:      sat - Nome convenzionale del satellite (e.g. L1)
% #              ver - Versione del modello
% #              t_camp - Periodo di campionamento (s)
% #              n_camp - Numero di campioni nel file
% #              t_inizio - Tempo di inizio nel formato ISO (e.g. 2001-01-01T00:00:00.000)
% #              t_fine - Tempo di fine nel formato ISO (e.g. 2001-01-01T00:00:00.000)
% #
% #              tempo_iso - Tempo del campione nel formato ISO (e.g. 2001-01-01T00:00:00.000)
% #              acc_x, acc_y, acc_z - Componenti del vettore accelerazione (m s^-2)
%
% # Metadata
% sat = <sat>
% ver = <version>
% t_camp = <periodo_campionamento>
% n_camp = <numero_campioni>
% t_inizio = <tempo_inizio_iso>
% t_fine = <tempo_fine_iso>
%
% # Data
% # t acc_x acc_y acc_z
% <tempo_iso_1> <acc_x_1> <acc_y_1> <acc_z_1>
% <tempo_iso_2> <acc_x_2> <acc_y_2> <acc_z_2>
% ...
% <tempo_iso_n> <acc_x_n> <acc_y_n> <acc_z_n>

path_save = 'D:\OneDrive - INAF - Istituto Nazionale di Astrofisica\G4S\Matlab_ray\Results\';
File_n = [desc{1,1},'_',nomefile,'.acc_est.asc'];
Filename = [path_save,File_n];

% Example description
description = {...
    ['# Filename:    ',File_n]...
    ['# Filetype:    ASCII'],...
    '# Author:      Carlo Lefevre (email:carlo.lefevre@inaf.it)',...
    ['# Last update: ',datestr(now,'ddd dd mmm yyyy HH.MM.SS')],...
    ['# Description: Geodyn external accelerations, ',Note],...
    [''],...
    '# Metadata',...
    };

header = {...
    [strcat("sat = ",desc{1,1})],...
    ['ver = ','001'],...
    [strcat("t_camp   = ",' ',string(desc{1,4}))],...
    [strcat("n_camp   = ",' ',string(length(t)))]...
    [strcat("t_inizio = ",' ',string(t(1)))],...
    [strcat("t_fine   = ",' ',string(t(end)))],...
    };


description2 = {['']...,
    '# Data',...
    '# t acc_x acc_y acc_z',...
    };

% Example data with four columns
data = {t,...
    SRP_ECI(1,:)',...
    SRP_ECI(2,:)',...
    SRP_ECI(3,:)'',...
    };



% Specify the file name
fileID = fopen(Filename, 'w');

% Write description lines
for i = 1:length(description)
    fprintf(fileID, '%s\n', description{i});
end

for i = 1:length(header)
    fprintf(fileID, '%s\n', header{i});
end

for i = 1:length(description2)
    fprintf(fileID, '%s\n', description2{i});
end

% Write data with four columns
for i = 1:size(data{1,1},1)
    fprintf(fileID, '%s %d %d %d\n', data{1, 1}(i), data{1, 2}(i), data{1, 3}(i), data{1, 4}(i));
end

% Close the file
fclose(fileID);


end

function save_data_geodine_quatrnion(desc,t,ACC,nomefile,Note)

% # Filename:    <sat>_<modello_acc>.att_est.asc
% # Filetype:    ASCII
% # Author:      <nome_autore_file> (email:<email_autore_file>)
% # Last update: <data_aggiornamento> # Preferibilmente nel formato 'ddd dd mmm yyyy, hh.mm.ss, CET' dato dal comando Unix 'date'
% #
% # Description: Geodyn external attitude
% #              <descrizione_contenuto_file>
% #              ...
% #
% # Format:      sat - Nome convenzionale del satellite (e.g. L1)
% #              ver - Versione del modello
% #              t_camp - Periodo di campionamento (s)
% #              n_camp - Numero di campioni nel file
% #              t_inizio - Tempo di inizio nel formato ISO (e.g. 2001-01-01T00:00:00.000)
% #              t_fine - Tempo di fine nel formato ISO (e.g. 2001-01-01T00:00:00.000)
% #
% #              tempo_iso - Tempo del campione nel formato ISO (e.g. 2001-01-01T00:00:00.000)
% #              q_1, q_2, q_3, q_4 - Componenti del quaternione assetto

path_save = 'C:\Users\Carlo\OneDrive - INAF - Istituto Nazionale di Astrofisica\G4S\Matlab_ray\Results\';
File_n = [desc{1,1},'_',nomefile,'.att_est.asc'];
Filename = [path_save,File_n];

% Example description
description = {...
    ['# Filename:    ',File_n]...
    ['# Filetype:    ASCII'],...
    '# Author:      Carlo Lefevre (email:carlo.lefevre@inaf.it)',...
    ['# Last update: ',datestr(now,'ddd dd mmm yyyy HH.MM.SS')],...
    ['# Description: Geodyn external accelerations, ',Note],...
    [''],...
    '# Metadata',...
    };

header = {...
    [strcat("sat = ",desc{1,1})],...
    ['ver = ','001'],...
    [strcat("t_camp   = ",' ',string(desc{1,4}))],...
    [strcat("n_camp   = ",' ',string(length(t)))]...
    [strcat("t_inizio = ",' ',string(t(1)))],...
    [strcat("t_fine   = ",' ',string(t(end)))],...
    };


description2 = {['']...,
    '# Data',...
    '# t q_1 q_2 q_3 q_4',...
    };

% Example data with four columns
data = {t(1:length(ACC.Nx(1,:)))',...
      zeros(length(ACC.Nx(1,:)),1),...
      zeros(length(ACC.Nx(1,:)),1),...
      zeros(length(ACC.Nx(1,:)),1),...
       ones(length(ACC.Nx(1,:)),1)...
    };



% Specify the file name
fileID = fopen(Filename, 'w');

% Write description lines
for i = 1:length(description)
    fprintf(fileID, '%s\n', description{i});
end

for i = 1:length(header)
    fprintf(fileID, '%s\n', header{i});
end

for i = 1:length(description2)
    fprintf(fileID, '%s\n', description2{i});
end

% Write data with four columns
for i = 1:size(data{1,1},2)
    fprintf(fileID, '%s %d %d %d %d\n', data{1, 1}(i), data{1, 2}(i), data{1, 3}(i), data{1, 4}(i), data{1, 5}(i));
end

% Close the file
fclose(fileID);


end

function utcTime = unixTimeToUTC(unixTime)
% Convert Unix time (seconds since epoch) to UTC time
utcTime = datetime(unixTime, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
end
