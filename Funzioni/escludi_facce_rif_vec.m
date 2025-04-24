function Dentro_1 = escludi_facce_rif_vec(sat0,Rest,Centri_mesh,N_ray_mesh,P_source,ray_dir,ray_face,varargin);
Dentro_1 = [];
tt = 1;

IDX_Facce_ill = Rest.facce_restanti;       % variabile dinamica, indici facce restanti
raggio_dir    = ray_dir;
N_dir         = N_ray_mesh;

Dentro(IDX_Facce_ill) = false;

idx_Face   = sat0.Faces(ray_face,:);         % indici riga vertici faccia del raggio sorgente

% posizioni vertici della sorgente in esame
A        = sat0.Vertex(idx_Face(1),:);
B        = sat0.Vertex(idx_Face(2),:);
C        = sat0.Vertex(idx_Face(3),:);

% punti mesh del piano da analizzare
C1_mesh       = sat0.Vertex(sat0.Faces(IDX_Facce_ill,1),:);       % Punto 1 mesh del piano illuminato
C2_mesh       = sat0.Vertex(sat0.Faces(IDX_Facce_ill,2),:);
C3_mesh       = sat0.Vertex(sat0.Faces(IDX_Facce_ill,3),:);

% proiezione della mesh del raggio sul piano in esame
N_mesh = sat0.Normals_mesh(IDX_Facce_ill,:);

rag_dir_n = repmat(raggio_dir,size(C1_mesh,1), 1);
% Calcola le intersezioni tra il raggio e il piano per tutte le facce rimanenti
A_prj2 = A + sum((rag_dir_n.*(A-C1_mesh)),2)./sum(rag_dir_n.*N_mesh,2).*rag_dir_n; % line_plane_intersection(raggio_dir, A, N_mesh, C1_mesh);
B_prj2 = B + sum((rag_dir_n.*(B-C2_mesh)),2)./sum(rag_dir_n.*N_mesh,2).*rag_dir_n;
C_prj2 = C + sum((rag_dir_n.*(C-C3_mesh)),2)./sum(rag_dir_n.*N_mesh,2).*rag_dir_n;

% Controlla se il raggio colpisce ciascuna faccia in esame
Dentro1 = any(ver_inner_vect(A_prj2, B_prj2, C_prj2, C1_mesh, cell2mat(varargin{1, 1})), 2);
Dentro2 = any(ver_inner_vect(A_prj2, B_prj2, C_prj2, C2_mesh, cell2mat(varargin{1, 1})), 2);
Dentro3 = any(ver_inner_vect(A_prj2, B_prj2, C_prj2, C3_mesh, cell2mat(varargin{1, 1})), 2);

% Trova gli indici delle facce che soddisfano le condizioni
if max(Dentro1) == 1
    Dentro_1 = IDX_Facce_ill(Dentro1);
elseif max(Dentro2) == 1
    Dentro_1 = IDX_Facce_ill(Dentro2);
elseif max(Dentro3) == 1
    Dentro_1 = IDX_Facce_ill(Dentro3);
end



end