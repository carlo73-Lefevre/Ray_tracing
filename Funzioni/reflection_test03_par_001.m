
function ray_out = reflection_test03_par_001(ray, sat0, tolerance)
    % Preallocazione dei risultati
    num_rays = length(ray.face_source);
    p_source_cell = cell(num_rays, 1);
    face_source_cell = cell(num_rays, 1);
    dir_cell = cell(num_rays, 1);
    face_from_cell = cell(num_rays, 1);

    % Estrazione dei dati necessari per ogni raggio
    P_rif = ray.p_source;
    rif = ray.dir;
    facce_source_idx = ray.face_source;

    parfor ii = 1:num_rays
        try
            % Calcolo degli indici delle facce illuminate
            idx_ill = indici_facce_ill(rif(ii,:), sat0, P_rif(ii,:), facce_source_idx(ii));

            if ~isempty(idx_ill)
                % Calcolo delle distanze e ordinamento delle facce potenziali
                VV_rest_sort = Calcola_distanze(sat0, idx_ill, facce_source_idx(ii));
                idx_facce_pot = VV_rest_sort(:,2);

                % Ricerca della faccia colpita
                [trovato, idx_riga_faccia, P_imp] = trova_faccia_colpita(sat0, idx_facce_pot, rif(ii,:), P_rif(ii,:), tolerance);

                if trovato
                    % Analisi delle facce complanari
                    [p_source_local, face_source_local, dir_local, face_from_local] = ...
                        analizza_facce_complanari(sat0, idx_riga_faccia, P_imp, rif(ii,:), facce_source_idx(ii));
                    
                    % Salvataggio dei risultati locali
                    p_source_cell{ii} = p_source_local;
                    face_source_cell{ii} = face_source_local;
                    dir_cell{ii} = dir_local;
                    face_from_cell{ii} = face_from_local;
                end
            end
        catch ME
            warning('Error in iteration %d: %s', ii, ME.message);
        end
    end

    % Combinazione dei risultati
    ray_out = c_Rays();
    ray_out.face_source = vertcat(face_source_cell{:});
    ray_out.p_source = vertcat(p_source_cell{:});
    ray_out.dir = vertcat(dir_cell{:});
    ray_out.face_from = vertcat(face_from_cell{:});
end

function [trovato, idx_riga_faccia, P_imp] = trova_faccia_colpita(sat0, idx_facce_pot, rif, P_rif, tolerance)
    trovato = false;
    idx_riga_faccia = [];
    P_imp = [];

    for ss = 1:length(idx_facce_pot)
        idx_riga_faccia = idx_facce_pot(ss);
        if idx_riga_faccia > size(sat0.Faces, 1)
            % warning('Invalid face index: %d', idx_riga_faccia);
            continue;
        end
        idx_righe_Vertex = sat0.Faces(idx_riga_faccia,:);
        A = sat0.Vertex(idx_righe_Vertex(1),:);
        B = sat0.Vertex(idx_righe_Vertex(2),:);
        C = sat0.Vertex(idx_righe_Vertex(3),:);

        if abs(dot(rif, sat0.Normals_mesh(idx_riga_faccia,:))) < tolerance
            continue;
        end

        P_imp = Intersection_CGPT(rif, P_rif, sat0.Normals_mesh(idx_riga_faccia,:), A);

        if norm(P_imp) < 1e5 && in_or_out(A, B, C, P_imp)
            trovato = true;
            break;
        end
    end
end

function [p_source, face_source, dir, face_from] = analizza_facce_complanari(sat0, idx_riga_faccia, P_imp, rif, facce_source_idx)
    if idx_riga_faccia > length(sat0.Face_coplanar)
        % warning('Invalid idx_riga_faccia: %d', idx_riga_faccia);
        p_source = []; face_source = []; dir = []; face_from = [];
        return;
    end

    Compl_facce_idx = sat0.Face_coplanar == sat0.Face_coplanar(idx_riga_faccia);
    idx_compl = find(Compl_facce_idx == 1);

    if isempty(idx_compl)
        % warning('No coplanar faces found for idx_riga_faccia: %d', idx_riga_faccia);
        p_source = []; face_source = []; dir = []; face_from = [];
        return;
    end

    Rest.facce_restanti = idx_compl;

    N = sat0.Normals_mesh(idx_riga_faccia,:);
    ray_rif_dir = 2 * (dot(N, -rif)) .* N + rif;
    
    idx_face_rif = escludi_facce_rif(sat0, Rest, sat0.Centers_mesh, N, P_imp, rif, facce_source_idx);
    
    if isempty(idx_face_rif)
        idx_face_rif = idx_riga_faccia;
    end

    ray_rif = repmat(ray_rif_dir, length(idx_face_rif), 1);
    P_imp = sat0.Centers_mesh(idx_face_rif,:);

    p_source = P_imp;
    face_source = idx_face_rif';
    dir = ray_rif;
    face_from = facce_source_idx * ones(size(P_imp,1),1);
end
