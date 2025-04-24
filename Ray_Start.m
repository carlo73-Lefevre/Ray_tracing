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

% path= 'C:\Users\carlo\OneDrive - INAF - Istituto Nazionale di Astrofisica\G4S\Matlab_ray';
path = 'D:\OneDrive - INAF - Istituto Nazionale di Astrofisica\G4S\Matlab_ray';

addpath(path);
addpath([path,'\Funzioni\'])
addpath([path,'\Dati Caomsol\'])
addpath([path,'\Class\'])

tic
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
Sat_inp.Box    = 'Box_WING_SAT';
Sat_inp.Wings  = 'Box_WING_Wings_center_003';

Sat  = SAT(); % Classe dati satellite statici
Sat1 = SAT(); % Classe dati satellite statici

% Parametri satellite ordinati per material idx
Sat1  = Read_mesh_MAT_txt(Sat_inp.Box,path);

if  length(fields(Sat_inp)) == 2
    Sat2 = SAT(); % Classe dati satellite statici
    Sat2  = Read_mesh_MAT_txt(Sat_inp.Wings,path);
end

% beta_y : 0 -> Z-   90 -> X-  (-90->X+ da mantenere all'ombra!!!)

% beta_y  = [linspace(0,10,20) ...
%     linspace(11,88,60) ...
%     linspace(88.1,90,10) ...
%     linspace(90.1,92,10) ...
%     linspace(93,170,30)...
%     linspace(171,180,20)]; % angolo Z intorno a Y
% 
% alf_x   = [linspace(90,100,20)...
%     linspace(100,178,60)...
%     linspace(178.1,180,10) ...
%     linspace(180,178.1,10)...
%     linspace(178.1,100,30)...
%     linspace(100,90,20)...
%     ]; %zeros(length(beta_y));           % angolo Z intorno a X

% [X_bet,Y_alf] = meshgrid(beta_y, alf_x);
% Z = (X_bet+Y_alf)*0;
% surf(X_bet,Y_alf,Z)

sun0_dir = [0 0 -1]; % [-1 -1 ];

Sat   = SAT();
ray   = c_Rays();

% ray_b(1:length(beta_y),1:5,length(alf_x))  = ray;
% Sat_b(1:length(beta_y),length(alf_x))      = Sat;

video = 'off';

if strcmp(video,'on')
    path_dave = ('D:\OneDrive - INAF - Istituto Nazionale di Astrofisica\G4S\Matlab_ray\');
    v = VideoWriter([path_dave,'G4s_BOX_W_vero.mp4'],'MPEG-4');
    open(v);
end

%% posizione del vettore sole da Calcolo_rad_sol_alb_IR_12x24_f_Carlo001
%  t_launch = utc2mjd([2015 12 17 11 51 00]);  
load('Sun_Body_dir.mat');

Xt = [ones(1,size(sun_body,2)).*[1 0 0]'];
Yt = [ones(1,size(sun_body,2)).*[0 1 0]'];
Zt = [ones(1,size(sun_body,2)).*[0 0 1]'];
 % angoli del vettore Sole nel ref Body
beta_sun = rad2deg(acos((dot(sun_body,Zt))));
alfa_sun = rad2deg(acos((dot(SUN,Xt))));
%%

tic
ACC.beta = -10; kk = 0;
for ii = 1:24*7*60 %length(beta_sun)

    if beta_sun(ii) < 12 || beta_sun(ii) > 168
         if min(abs(ACC.beta-beta_sun(ii))) < 0.005 %risoluzione angolare
            continue
         end
    elseif beta_sun(ii) > 88 || beta_sun(ii) < 92
         if min(abs(ACC.beta-beta_sun(ii))) < 0.005 %risoluzione angolare
            continue
         end
    else 
         if min(abs(ACC.beta-beta_sun(ii))) < 0.3 %risoluzione angolare
            continue
         end
    end
    % unisce i 3D ruotando i pannelli e calcola le facce complanari
    if exist('Sat2')
        [Sat,Ry] = monta_sat(Sat1,Sat2,-beta_sun(ii));
    else
        Sat = Sat1;
        Sat = trova_complanari(Sat);
        Sat = centri_mesh(Sat);

        Ry = [cos(deg2rad(-beta_sun(ii)))  0  sin(deg2rad(-beta_sun(ii)));...
                 0                         1          0               ;...
            -sin(deg2rad(-beta_sun(ii)))   0 cos(deg2rad(-beta_sun(ii)))];
    end


   kk = kk+1;

    Rx = [1                0                0                  ; ...
        0    cos(deg2rad(alfa_sun(kk))) -sin(deg2rad(alfa_sun(kk)));...
        0    sin(deg2rad(alfa_sun(kk)))  cos(deg2rad(alfa_sun(kk)))];

    % direzione del Sole
    ray_dir_a     =  [Rx*Ry*sun0_dir']';
    sun_dir       =  -SUN(:,ii)';% ray_dir_a ./ norm(ray_dir_a);  % coseni direttori raggi incidenti


    [Sat,ray,frame]    = Ray_tracing(Sat,sun_dir,video);

    if strcmp(video,'on')
        writeVideo(v,frame);
    end

    [ACC_S1,ACC_N1,Area] = SRP(Sat,ray(1,1),ray(1,2),ray(1,1).dir(1,:));
   
    ACC.beta(kk)        = beta_sun(ii);
    ACC.alf_x(kk)       = alfa_sun(ii);
    ACC.sun(kk,:)       = ACC_S1;
    ACC.sun_mag(kk)     = norm(ACC_S1);

    ACC.N(kk,:)         = ACC_N1;
    ACC.N_mag(kk,:)     = norm(ACC_N1);

    ACC.SRP_mag(kk)      = norm(ACC_S1 + ACC_N1 );

    Sat_b(kk) = Sat;

    % popola ray_b tenendo conto che ray può andara da 1 a 5
    for jj =1:size(ray,2)
        ray_b(kk,jj) = ray(jj);
    end

end



toc

if strcmp(video,'on')
    close(v);
end

toc

plot(ACC.beta,ACC.SRP_mag,'.-')


%%
% load('Sun_Body.mat');
% Xt = [ones(1,length(SUN)).*[1 0 0]'];
% Yt = [ones(1,length(SUN)).*[0 1 0]'];
% Zt = [ones(1,length(SUN)).*[0 0 1]'];
%
% beta_sun = rad2deg(acos((dot(SUN,Zt))));
% alfa_sun = rad2deg(acos((dot(SUN,Xt))));


ACC_beta_abs_time = interp1(ACC.beta,ACC.SRP_mag,beta_sun,'spline');
ACC_alfa_abs_time = interp1(ACC.alf_x,ACC.SRP_mag,alfa_sun,'spline');
% ACC_abs_time = interp2(ACC.beta,ACC.alf_x,ACC.SRP_mag,beta_sun,alfa_sun,'spline');

plot(t,ACC_beta_abs_time)
plot(t,ACC_alfa_abs_time)
%%

%  grafica_riflessioni(Sat_b(1),ray_b(1,2))
%  grafica_riflessioni(Sat_b(1),ray_b(1,3))
%  grafica_riflessioni(Sat_b(1),ray_b(1,4))
%  grafica_riflessioni(Sat_b(1),ray_b(1,5))

X = [1 0 0];Y =[0 1 0]; Z=[0 0 1];

figure(1)
clf;
grafica_riflessioni(Sat_b(1),ray_b(1,1))
quiver3([0 0 0],[0 0 0],[0 0 0],X,Y,Z,6);
view(130,30)
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
% Sat_materiali(Sat)

%
figure(2)
clf
grafica_riflessioni(Sat,ray_b(1,1))

figure(3)
clf;
grafica_riflessioni(Sat_b(1),ray_b(1,2));


figure(5)
clf
grafica_riflessioni_source(Sat,ray_b(1,2),ray_b(1,3))

figure(6)
clf
grafica_riflessioni_source(Sat,ray_b(1,3),ray_b(1,4))

figure(7)
clf
grafica_riflessioni_source(Sat,ray_b(1,4),ray_b(1,5))
%%
Sat_materiali(Sat)

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





%% colora facce materiali
% %
function Sat_materiali(Sat)
figure(2)
clf
hold on;
axis equal
box on;

p2 = patch('Faces',Sat.Faces,'Vertices',Sat.Vertex,'EdgeColor','none'); % facce colpite dai riflessi
p2.FaceColor = [0.7,0.7,0.7];
hold on

cmap =[0.5 0.5 0.5;...
    0 1 0;...
    0 0 1;...
    1 0 0;...
    0.0 0.5 0.5;...
    0.5 0.5 0.0;...
    0.5 0.5 1;...
    0.1 0.5 1;];

L ={};
for ii = 0:max(Sat.materiali_idx) +1
    idx_mat_log = (Sat.materiali_idx == ii);  % vertici che contengono il mat
    [idx_row,idx_col] = find(idx_mat_log == 1);

    idx_row_Fn = [];
    for jj=1:length(idx_row)

        [idx_row_F,idx_col_F] = find(Sat.Faces == idx_row(jj));
        idx_row_Fn = [idx_row_Fn;idx_row_F];


    end

    p5 = patch('Faces',Sat.Faces(idx_row_Fn,:),'Vertices',Sat.Vertex,'EdgeColor','none'); % faccia in esame per capire se sta dietro
    p5.FaceColor = cmap(ii+1,:);
    p5.EdgeLighting =  'gouraud';
    axis equal

    h1(ii+1) = p5 ;

    L{ii+1,1} = {(Sat.optical_par(ii+1).name)};

end
legend(h1,string(L))
title('MATLAB SATELLITE Materials and Mesh','Interpreter','latex')

end

function [Sat,Ry] = monta_sat(Sat1,Sat2,beta_y)
% funzione che unisce i 2 3D, ruota i pannelli solari di un angolo beta_y
% intorno a Y
% e calcola quale faccesono complanari

% ruota pannelli di un angolo alfa
angy = deg2rad(beta_y);
% angz = deg2rad(0);

% Rx = [1         0        0; ...
%     0 cos(angx) -sin(angx);...
%     0 sin(angx) cos(angx)];

Ry = [cos(angy)  0  sin(angy);...
    0          1          0;...
    -sin(angy) 0 cos(angy)];

% Rz = [cos(angz)  -sin(angz)   0;...
%     sin(angz)   cos(angz)   0;...
%     0           0           1];


% ruota i pannelli solari
Sat2.Vertex = [Ry*Sat2.Vertex']';

% unisci i due body
Sat = unisci_bodys(Sat1,Sat2);

% determina i centri Mesh (da fare dopo la rotazione=)
% Sat = centri_mesh(Sat);

% determina quali facce sono complanari
% L'induce riga è relativo all'indice della Facce in esame che sono
% complanari se il valore è uguale
% idx_faces = [false(length(Sat.Faces),1)]';
% App = (1:length(Sat.Faces))';
% App = cerca_piani_app(idx_faces,Sat.Centri_mesh,App,1,Sat);
% Sat.Face_rige_compl = App;

Sat = trova_complanari(Sat);

end

function Sat = trova_complanari(Sat)

% determina quali facce sono complanari
% L'induce riga è relativo all'indice della Facce in esame che sono
% complanari se il valore è uguale

%  trova i centri delle mesh
Sat       = centri_mesh(Sat);

idx_faces = [false(length(Sat.Faces),1)]';
App       = (1:length(Sat.Faces))';

% cerca le mesh complanari
App = cerca_piani_app(idx_faces,Sat.Centri_mesh,App,1,Sat);

Sat.Face_rige_compl = App;

end
