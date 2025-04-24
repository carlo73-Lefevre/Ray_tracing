function [face] = Face_illuminated(MAT_tot,ray_dir,ray_start)

% Funzione che ciclando per tutte le facce della mesh determina se la faccia è colpita
% dal raggio escludendo tutte quelle con una normale maggiore di 90°
% (angolo conico) rispetto al raggio

face = [];
%% Triangolo relativo alla facce F2(1,:)
% Dentro = false;

for Zona = 1 : size(MAT_tot,1) % cicla per facce con stesso materiale
    MAT = MAT_tot(Zona);       % inizializza la MAT della zona
    for Fnum = 1:size(MAT.Face,1) % cicla per tutte le facce

        % Nodi della faccia in esame
        A = MAT.Vertex(MAT.Face(Fnum,1),:);
        B = MAT.Vertex(MAT.Face(Fnum,2),:);
        C = MAT.Vertex(MAT.Face(Fnum,3),:);

        %% calcola i punto di impatto
        v1 = vector(A,B);
        v2 = vector(B,C);
        v3 = vector(C,A);

        Plane_norm = cross(v1,v3);   % Normal to plane
        M = A;                       % Point in the plane

        if rad2deg(acos(dot(Plane_norm/norm(Plane_norm),ray_dir/norm(ray_dir)))) < 90 % verifica che il raggio arrivi davanti

            [P_ray_imp,~] = line_plane_intersection(ray_dir, ray_start, Plane_norm, M); % intersezione Piano Linea

            % Verifica che il punto sia all'interno
            if isempty(P_ray_imp)  % se il raggio è complanare P_ray_imp = []

                Dentro = false;  

            else

                Dentro = ver_inner(v1,v2,v3,A,B,C,P_ray_imp);
                
                %%  verifica che la normale al piano sia opposta
                if Dentro == true
                    n          = Plane_norm/norm(Plane_norm);                                % narmale al piano di impatto
                    d          = dot(ray_dir/norm(ray_dir),Plane_norm/norm(Plane_norm));     % angolo tra piano e raggio
                    ray_rif_dir = ray_dir -2*(dot(ray_dir,n)).*n;

                    path_ray = norm(ray_start -P_ray_imp);                                   % distanza P impatto - punto sorgente


                    face_n.area               = calc_area(v1,v2,v3);
                    face_n.Material_index     = MAT.idx;
                    face_n.vertex_A           = A;
                    face_n.vertex_B           = B;
                    face_n.vertex_C           = C;
                    face_n.view_fact          = d; % cos beta
                    face_n.impact_P           = P_ray_imp;
                    face_n.dist               = v_leng(P_ray_imp - ray_start);
                    face_n.ray_source_P_leng  = round(v_leng(ray_start),5);
                    face_n.ray_source_P       = ray_start;
                    face_n.ray_source_dir     = ray_dir;
                    face_n.ray_rif_dir        = ray_rif_dir;
                    face_n.Face               = MAT.Face(Fnum,:);

                    face      = [face face_n];

                end
            end
        end
    end
end

end






