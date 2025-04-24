%%                      Ray Tracing Comsol - Matlab
%                    Coded by Carlo Lefevre 27-09-2022
% Il prog prende come input i file esportati da Comsol. Il .dat è un testo che
% traccia i parametri ottici, il .xml identifica quali materiali uano i
% parametri ottici definiti in dat., nel .txt ci sono le matrici posizione
% dei nodi della mesh, la matrice posizioni delle facce di mesh (ogni riga
% ha 3 numeri che sono gli indici di riga dela matrice dei nodi) e un
% ultimo vettore che associa ad ogni facci quale materiale è stato usato
% (con un indice che va da 0 a n-1 materiali).

function [Sat,ray,frame] = Ray_tracing(Sat,sun_dir,ray1_rif,video)
% classe Raggi
ray0          = c_Rays();

n_dir =  repmat(sun_dir,size(Sat.Centers_mesh,1), 1);     % direzione dei raggi
CM    =  Sat.Centers_mesh;                                % centri delle mesh
u_dir =  repmat(-sun_dir,size(Sat.Centers_mesh,1), 1);    % normale al piano dove proiettare i centri mesh
PN    =  repmat(-sun_dir*5,size(Sat.Centers_mesh,1), 1);  % punti appartenenti al piano dove si vuole proiettare (sempre lo stesso)

% t = n*(M-N) / (u*n) formula che clcolala distanza lungo N-M  tra 
t = sum(n_dir.*(PN-CM),2) ./ sum(n_dir .* u_dir,2);

% intersection=N+t⋅u  formula che trova il punto di intersezione sul piano 
ray0.p_source = CM + t.*u_dir;

% direzioni raggio dal piano sorgente
ray0.dir  = -u_dir;

% Trova le normali alla mesh e i raggi riflessi (tutti)
ray1 = c_Rays();
% [Sat.Normals_mesh ,ray1_dir] = trova_Normali(Sat,ray0.dir);

% trova le facce potenzialmente illuminate ordinate dalla più vicina 
[MAT,loc]         =  facce_illuminate_potenziali(Sat,ray0);
ray0.p_source     =  ray0.p_source(MAT.ordine_per_distanza,:);
ray0.dir          =  ray0.dir(MAT.ordine_per_distanza,:);

% grafica facce restanti
% figure(21)
% clf;
% p0 = patch('Faces',Sat.Faces,'Vertices',Sat.Vertex); % facce colpite dai riflessi
% p0.FaceColor = 'w';
% axis equal;
% hold on
% p1 = patch('Faces',Sat.Faces(MAT.Dis_sort_ill,:),'Vertices',Sat.Vertex); % facce colpite dai riflessi
% p1.FaceColor = [0.2,0.2,1];
% view(-25,-20);

%%% *Funzione ricorsiva di eliminazione delle facce in ombra %%%

 MAT.Dis_sort_ill = escludi_facce(Sat,MAT,ray0,1);       % indici delle facce restanti

% grafica facce restanti
% figure(22)
% clf;
% p0 = patch('Faces',Sat.Faces,'Vertices',Sat.Vertex); % facce colpite dai riflessi
% p0.FaceColor = [0.7,0.7,0.7];
% axis equal;
% hold on
% p1 = patch('Faces',Sat.Faces(MAT.Dis_sort_ill,:),'Vertices',Sat.Vertex); % facce colpite dai riflessi
% p1.FaceColor = [0.2,0.2,0.55];


MAT.source         = MAT.Dis_sort_ill ;                 % facce illuminate restanti;
MAT.n_riflesso     = ones(length(MAT.Dis_sort_ill),1);  % variabile che tiene traccia del numero di riflessioni

% ricava facce in ombra
MAT.facce_scure  = 1 : size(Sat.Faces,1) ;
MAT.facce_scure(MAT.facce_restanti) = []; % tutte le facce - le illuminate

%%%%      %%% Raggi riflessi %%%     %%%%

% prima riflessione

ray1.p_source    = Sat.Centers_mesh(MAT.source,:);    % centri facce raggi riflessi
ray1.face_source = MAT.source;                       % indici delle facce dei riflessi
ray1.dir         = ray1_rif(MAT.source,:);           % direzione raggi riflessi

if size(ray1.p_source,1) > 0
frame = 0;
    if strcmp(video,'on')
        f = figure('visible','off');
        clf;
        grafica_riflessioni(Sat,ray1);
        title('Ray 1')
        xlim([-3/6 3/6])
        ylim([-7/14 7/14])
        zlim([-1/3 1/3])
        view(9,+32)

        quiver3(0,0,0,sun_dir(end,1),sun_dir(end,2),sun_dir(end,3),8); % direzione sole

        axis tight manual
        set(gca,'nextplot','replacechildren');
        frame = getframe(gcf);
    end


    %% Seconda riflessione
    
      % ray2 = reflection_test03(ray1,Sat,1e-7^2); %
     
    % % 
    % if size(ray2.face_source,1) > 0
    % 
    %% Terza riflessione
    % 
    %     ray3 = reflection_test03(ray2,Sat,1e-6);
    % 
    % 
         %% Quarta riflessione
    % 
    %     if size(ray3.p_source,1) > 0
    % 
    %         ray4  = reflection_test03(ray3,Sat,1e-6);
    % 
    %     end
    % end
end

ray = ray0;

if exist('ray1','var')
    ray = [ray,ray1];

    if exist('ray2','var') && ~isempty(ray2.face_source)
        ray = [ray,ray2];

        if exist('ray3','var') && ~isempty(ray3.face_source)
            ray = [ray,ray3];

            if exist('ray4','var') && ~isempty(ray4.face_source)
                ray = [ray,ray4];
            end
        end
    end
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
% %  plot_quiv(ray1.p_source,ray1.dir) % raggi riflessi che colpiscono una faccia
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




end