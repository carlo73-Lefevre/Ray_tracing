%%                      Ray Tracing Comsol - Matlab
%                    Coded by Carlo Lefevre 27-09-2022
% Il prog prende come input i file esportati da Comsol. Il .dat è un testo che
% traccia i parametri ottici, il .xml identifica quali materiali uano i
% parametri ottici definiti in dat., nel .txt ci sono le matrici posizione
% dei nodi della mesh, la matrice posizioni delle facce di mesh (ogni riga
% ha 3 numeri che sono gli indici di riga dela matrice dei nodi) e un
% ultimo vettore che associa ad ogni facci quale materiale è stato usato
% (con un indice che va da 0 a n-1 materiali).

function [Sat,Ray_tracing(beta)
% close all
% path= 'C:\Users\carlo\OneDrive - INAF - Istituto Nazionale di Astrofisica\G4S\Matlab_ray';
path = 'D:\OneDrive - INAF - Istituto Nazionale di Astrofisica\G4S\Matlab_ray';
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
% retro_rif retro_rif_low  [1 -1 -1];
% bullone
% specchio  [0 1 0]
% Cannon Cannon_fine [0 1 -1];
% Lageos_Carlo Lageos_Carlo_Mas Lageos_Carlo_Mas2
% Box_WING_SAT Box_WING_Wings
% 
% Sat  = SAT(); % Classe dati satellite statici
% Sat1 = SAT(); % Classe dati satellite statici
% Sat2 = SAT(); % Classe dati satellite statici

% Parametri satellite ordinati per material idx
% Sat1  = Read_mesh_MAT_txt('Box_WING_SAT',path);
% Sat2  = Read_mesh_MAT_txt('Box_WING_Wings',path);


% ruota pannelli di un angolo alfa
angx = deg2rad(0);
angy = deg2rad(0);
angz = deg2rad(0);

Rx = [1         0        0; ...
    0 cos(angx) -sin(angx);...
    0 sin(angx) cos(angx)];

Ry = [cos(angy)  0  sin(angy);...
    0          1          0;...
    -sin(angy) 0 cos(angy)];

Rz = [cos(angz)  -sin(angz)   0;...
    sin(angz)   cos(angz)   0;...
    0           0           1];

% ruota i pannelli solari
Sat2.Vertex = [Ry*Sat2.Vertex']';

% unisci i due body
Sat = unisci_bodys(Sat1,Sat2);

% direzione del Sole
ray_dir_a     =  [Rx*Ry*Rz*[0 1 0]']';
sun_dir       = ray_dir_a ./ norm(ray_dir_a);  % coseni direttori raggi incidenti

% determina i centri Mesh (da fare dopo la rotazione=)
Sat = centri_mesh(Sat);

% classe Raggi
ray0          = c_Rays();

% Genera il Piano sorgente
for tt = 1 : size(Sat.Centri_mesh,1)
    % punti proiezione dei Vertici su un piano
    [ray0.p_source(tt,:),~] = line_plane_intersection(sun_dir, Sat.Centri_mesh(tt,:),sun_dir,  -sun_dir*8);

end


% direzioni raggio dal piano sorgente
ray0.dir  = sun_dir.*ones(size(ray0.p_source,1),3);

% Trova le normali alla mesh e i raggi riflessi (tutti)
ray1 = c_Rays();
[Sat.Normals_mesh ,ray1_dir] = trova_Normali(Sat,ray0.dir);


% trova le facce potenzialmente illuminate ordinate dalla più vicina
[MAT] = facce_illuminate_potenziali(Sat,ray0);

%%% *Funzione ricorsiva di eliminazione delle facce in ombra %%%

MAT.facce_restanti = escludi_facce(Sat,MAT,ray0,1);       % indici delle facce restanti

MAT.source         = MAT.facce_restanti ;                 % prime facce illuminate;
MAT.n_riflesso     = ones(length(MAT.facce_restanti),1);  % variabile che tiene traccia del numero di riflessioni

% ricava facce in ombra
MAT.facce_scure  = 1 : size(Sat.Faces,1) ;
MAT.facce_scure(MAT.facce_restanti) = []; % tutte le facce - le illuminate


             %%%%      %%% Raggi riflessi %%%     %%%%

% prima riflessione

ray1.p_source    = Sat.Centri_mesh(MAT.facce_restanti,:);    % centri facce raggi riflessi
ray1.face_source = MAT.facce_restanti;
ray1.dir         = ray1_dir(MAT.facce_restanti,:);
ì
grafica_riflessioni(Sat,ray1)
title('Ray 1')

% calcolo delle accelerazioni SRP dalla prima riflessione
[ACC_S,ACC_N] = SRP(Sat,ray0,ray1,sun_dir)

% Aeff = pi*0.1497^2;
% Cann_ball_err = 100*(norm(ACC_S) - (1380*Aeff/405.38)/299792458)/norm(ACC_S)


%% Seconda riflessione
% Da implementare la zona di sorgente come proiezione delle facce riflesse

if size(ray1.p_source,1) > 0

    % determina quali facce sono complanari
    % L'induce riga è relativo all'indice della Facce in esame che sono
    % complanari se il valore è uguale
    idx_faces = [false(length(Sat.Faces),1)]';
    App = (1:length(Sat.Faces))';
    App = cerca_piani_app(idx_faces,Sat.Centri_mesh,App,1,Sat);
    Sat.Face_rige_compl = App;

    ray1_u = c_Rays();

    [ray1_u.face_source idx_u1 ] =  unique(ray1.face_source);  % elimina le facce contate più volte
    ray1_u.p_source    =  ray1.p_source(idx_u1,:);
    ray1_u.dir         =  ray1.dir(idx_u1,:);

%     [MAT,ray2] = reflection_test03(MAT,ray1,Sat,0.001^2); %

    ray2_u = c_Rays();

    [ray2_u.face_source idx_u2 ] =  unique(ray2.face_source);  % elimina le facce contate più volte
    ray2_u.p_source    =  ray2.p_source(idx_u2,:);
    ray2_u.dir         =  ray2.dir(idx_u2,:);
    ray2_u.face_from   =  ray2.face_from(idx_u2,:);


    p1 = grafica_riflessioni(Sat,ray2);
    title('Ray 2');
    p1.FaceColor = 'c';

%     Link = linkprop([ax1, ax2],{'CameraUpVector', 'CameraPosition', 'CameraTarget'});
%     setappdata(gcf, 'StoreTheLink', Link);
%     view(angz,angy)

end

if size(ray2_u.face_source,1) > 0

    [ACC_S2,ACC_N2] = SRP(Sat,ray1_u,ray2_u,sun_dir)


    %% Terza riflessione

    if size(ray2_u.p_source,1) > 0

        [MAT,ray3] = reflection_test03(MAT,ray2,Sat,1e-6);

        ray3_u = c_Rays();

        [ray3_u.face_source idx_u3 ] =  unique(ray3.face_source);  % elimina le facce contate più volte
        ray3_u.p_source    =  ray3.p_source(idx_u3,:);
        ray3_u.dir         =  ray3.dir(idx_u3,:);
        ray3_u.face_from   =  ray3.face_from(idx_u3,:);

        p1 = grafica_riflessioni(Sat,ray3);
        title('Ray 3');
        p1.FaceColor = 'r';

    end

    if size(ray3_u.face_source,1) > 0

        %         [ACC_S3,ACC_N3] = SRP(Sat,ray2_u,ray3_u,sun_dir)

        %% Quarta riflessione

        if size(ray3_u.p_source,1) > 0

            [MAT,ray4] = reflection_test03(MAT,ray3,Sat,1e-6);

            ray4_u = c_Rays();

            [ray4_u.face_source idx_u4] =  unique(ray4.face_source);  % elimina le facce contate più volte
            ray4_u.p_source    =  ray4.p_source(idx_u4,:);
            ray4_u.dir         =  ray4.dir(idx_u4,:);
            ray4_u.face_from   =  ray4.face_from(idx_u4,:);

            p1 = grafica_riflessioni(Sat,ray4);
            title('Ray 4');
            p1.FaceColor = 'm';

            %             [ACC_S4,ACC_N4] = SRP(Sat,ray3_u,ray4_u,sun_dir)

        end
    end
end

ray = [ray0,ray1,ray2,ray4];


end

%%  %%% Zona grafici %%%






% da verificare come graficare solo i raggi riflessi ray1 che colpiscono
% le facce di ray2
%  figure(21)
%  clf;
%  hold on;
%  p0 = patch('Faces',Sat.Faces,'Vertices',Sat.Vertex); % facce colpite dai riflessi
%  p0.FaceColor = [0.7,0.7,0.7];
%  p01 = patch('Faces',Sat.Faces(ray2.face_source,:),'Vertices',Sat.Vertex); % facce colpite da ray1
%  p01.FaceColor = 'y';
% %  ray1_face     = ray1.face_source(ismember(ray1.face_source,ray2.face_source),:);
% %  dir_colp      = ray1.dir(ismember(ray1.face_source,ray1_face),:);
% %  p_source_colp = ray1.p_source(ismember(ray1.face_source,ray1_face),:);
%  axis equal;
%  box on;
%  plot_quiv(ray1.p_source,ray1.dir) % raggi riflessi che colpiscono una faccia
%



%  graph_riflessi(Sat,MAT,ray2)
% title('Ray 2')
%  p2 = patch('Faces',Sat.Faces(ray3.face_source,:),'Vertices',Sat.Vertex); % facce colpite dai riflessi
%  p2.FaceColor = 'r';
%
%
%  graph_riflessi(Sat,MAT,ray3)
%   p4 = patch('Faces',Sat.Faces(ray4.face_source,:),'Vertices',Sat.Vertex); % facce colpite dai riflessi
%  p4.FaceColor = 'r';
%   title('Ray 3')
%  p1 = patch('Faces',Sat.Faces(ray3.face_source,:),'Vertices',Sat.Vertex);
%  p1.FaceColor = 'r';



%% colora facce materiali

% 
% figure(7)
% clf
% hold on;
% axis equal
% box on;
% 
% p2 = patch('Faces',Sat.Faces,'Vertices',Sat.Vertex,'EdgeColor','none'); % facce colpite dai riflessi
% p2.FaceColor = [0.7,0.7,0.7];
% hold on
% for ii = 0:max(Sat.materiali_idx) +1
%     idx_mat_log = (Sat.materiali_idx == ii);  % vertici che contengono il mat
%     [idx_row,idx_col] = find(idx_mat_log == 1);
%     col(ii+1,:) =[rand rand rand];
%     for jj=1:length(idx_row)
% 
%         [idx_row_F,idx_col_F] = find(Sat.Faces == idx_row(jj));
%         p5 = patch('Faces',Sat.Faces(idx_row_F,:),'Vertices',Sat.Vertex,'EdgeColor','none'); % faccia in esame per capire se sta dietro
%         p5.FaceColor = col(ii+1,:);
%         p5.EdgeLighting =  'gouraud';
%     end
% end
% 
% title('MATLAB SATELLITE Materials and Mesh','Interpreter','latex')
