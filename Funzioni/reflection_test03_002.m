
%% Calcola raggio riflesso e indice faccia illuminata

function ray_out = reflection_test03 (ray,Sat,varargin)

P_rif            = ray.p_source;    % punti partenza raggio (indice riga MAT.Centres)
rif_dir          = ray.dir;         % direzione raggi

facce_source_idx = ray.face_source; % indice facce sorgenti  (sat0.Faces(idx))

ray_out  = c_Rays(); % crea la classe ray
c = 1;

for ii = 1:length(facce_source_idx)  % cicla per ogni raggio riflesso (delle facce ill)

    % indice riga delle facce la cui normale è < 0 con la normale alla
    % faccia mesh e faccia posta davanti

    idx_ill = indici_facce_ill(rif_dir(ii,:),Sat,P_rif(ii,:),facce_source_idx(ii));

    % ordina le mesh dalla più vicina
    if ~isempty(idx_ill)
        VV_rest_sort  = Calcola_distanze(Sat,idx_ill,facce_source_idx(ii));
        % indici delle facce potenziali ordinate per distanza
        idx_facce_pot = (VV_rest_sort(:,2));

    else
        continue
    end


    % Passi:
    % 1. proietto su di un piano parallelo alla mesh di partenza tutte le mesh
    % ordinate per distanza

    % posizioni vertici della mesh in esame (per verifiace se è in ombra)
    idx_P    = Sat.Faces(idx_facce_pot,:);

    % Punti della mesh sul satellite (tra quelli illuminabili)
    A        = Sat.Vertex(idx_P(:,1),:);
    B        = Sat.Vertex(idx_P(:,2),:);
    C        = Sat.Vertex(idx_P(:,3),:);

    F_A = proietta_piano_mesh(P_rif(ii,:),rif_dir(ii,:),A,idx_facce_pot);
    F_B = proietta_piano_mesh(P_rif(ii,:),rif_dir(ii,:),B,idx_facce_pot);
    F_C = proietta_piano_mesh(P_rif(ii,:),rif_dir(ii,:),C,idx_facce_pot);


    % 2. trovo le  mesh che ha un punto inteno
    % indice di idx_facce_pot che indica le facce colpite
    Dentro_idx_A =  dentro_list_vect_002(F_A,F_B,F_C,P_rif(ii,:),Sat);

    if ~isempty(Dentro_idx_A) && abs(dot(Sat.Normals_mesh(idx_facce_pot(Dentro_idx_A(1)),:),rif_dir(ii,:))) > 1e-6
         
        dir_ref = calcola_raggio_riflesso(rif_dir(ii,:),Sat.Normals_mesh(idx_facce_pot(Dentro_idx_A(1)),:));


        ray_out.face_source(c,:)     = idx_facce_pot(Dentro_idx_A(1));     % indice riga delle Face
        ray_out.p_source(c,:)        = Sat.Centers_mesh(idx_facce_pot(Dentro_idx_A(1)),:) ;
        ray_out.dir(c,:)             = dir_ref  ;

        ray_out.face_from(c,:)       = ray.face_source(ii);
        ray_out.p_from(c,:)          = Sat.Centers_mesh(ray.face_source(ii),:);
        ray_out.dir_from(c,:)        = rif_dir(ii,:)  ;

        c = c+1;
      
    end
end

end

function F_A = proietta_piano_mesh(P_rif,raggio_dir,A,idx_facce_pot)
% proietta A su un piano che passa per P_rif e ha giacitura raggio_dir

n_dir =  repmat(raggio_dir,size(idx_facce_pot,1), 1);

% t = n*(M-N) / (u*n) formula che clcolala distanza lungo N-M  tra
t_F_A = sum(n_dir.*(P_rif-A),2) ./ sum(n_dir .* n_dir,2);


% intersection=N+t⋅u  formula che trova il punto di intersezione sul piano
% sorgente della proiezione dei punti di mesh
% coordinate di tutti i punti del Sat proiettati sul piano mesh
F_A = A + t_F_A.*n_dir;

end

function raggio_riflesso = calcola_raggio_riflesso(n_dir, normal)
% Normalizza i vettori
n_dir = n_dir / norm(n_dir);
normal = normal / norm(normal);

% Calcola il prodotto scalare
dot_product = dot(n_dir, normal);

% Calcola il raggio riflesso
raggio_riflesso = n_dir - 2 * dot_product * normal;
end





