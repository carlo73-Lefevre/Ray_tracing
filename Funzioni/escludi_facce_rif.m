function Dentro_1 = escludi_facce_rif(sat0,Rest,Centri_mesh,N_ray_mesh,P_source,ray_dir,ray_face,varargin)
% Restituisce tutti i gli indici delle facce complanari alla colpita che si
% trovano allinterno della proiezione della faccia sorgende del raggio.
% varargin è un valore di soglia

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



for count = 1:length(IDX_Facce_ill)

    % punti mesh del piano da analizzare
    C1_mesh       = sat0.Vertex(sat0.Faces(IDX_Facce_ill(count),1),:);       % Punto 1 mesh del piano illuminato
    C2_mesh       = sat0.Vertex(sat0.Faces(IDX_Facce_ill(count),2),:);
    C3_mesh       = sat0.Vertex(sat0.Faces(IDX_Facce_ill(count),3),:);

    % proiezione della mesh del raggio sul piano in esame
    N_mesh = sat0.Normals_mesh(IDX_Facce_ill(count),:);

    [A_prj(1,:),~] = line_plane_intersection(raggio_dir,A,N_mesh,C1_mesh);
    [B_prj(1,:),~] = line_plane_intersection(raggio_dir,B,N_mesh,C2_mesh);
    [C_prj(1,:),~] = line_plane_intersection(raggio_dir,C,N_mesh,C3_mesh);


    % proiezioni del triangolino fem sul piano sorgente
    Faccia_idx   = IDX_Facce_ill(count);   % questo indice è lo steso pe le Facce per i ray_dir_v e ray_start_p

    % direzione normale  piano raggio
    C_mesh     = Centri_mesh(Faccia_idx,:);       % Centri mesh da controllare

    % verifica quali Punti mesh del piano in esame stanno nella proiezione
    % della mesh raggio sul piano

    % verifica se il centro mesh sta dentro
    Dentro(count) = ver_inner(A_prj,B_prj,C_prj,sat0.Centers_mesh(IDX_Facce_ill(count),:),cell2mat(varargin{1,1}));

    % verifica se uno dei vertici sta dentro
    if Dentro(count) == false

        Dentro(count) = ver_inner(A_prj,B_prj,C_prj,C1_mesh,cell2mat(varargin{1,1}));

    end

    if Dentro(count) == false

        Dentro(count) = ver_inner(A_prj,B_prj,C_prj,C2_mesh,cell2mat(varargin{1,1}));

    end

    if Dentro(count) == false

        Dentro(count) = ver_inner(A_prj,B_prj,C_prj,C3_mesh,cell2mat(varargin{1,1}));

    end



    if Dentro(count) == true
     

        Dentro_1(tt) =  Faccia_idx;


        % se il punto di sovraposizione corrisponde la punto di mesh prj
        % eliminalo
         mesh_dist =[abs(norm(A_prj - C1_mesh))<1e-12,abs(norm(A_prj - C2_mesh))<1e-12,abs(norm(A_prj - C3_mesh))<1e-12,...
                     abs(norm(B_prj - C1_mesh))<1e-12,abs(norm(B_prj - C2_mesh))<1e-12,abs(norm(B_prj - C3_mesh))<1e-12,...
                     abs(norm(C_prj - C1_mesh))<1e-12,abs(norm(C_prj - C2_mesh))<1e-12,abs(norm(C_prj - C3_mesh))<1e-12];

        if      any(mesh_dist > 0,'all')

                Dentro_1(tt)  = [];
                Dentro(count) = false;
                tt = tt-1;
        end

        tt = tt+1;

    end


    if sum(Dentro) < 1
        Dentro_1 = [];
        tt = 1;
    end




end