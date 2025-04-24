
%% Calcola raggio riflesso e indice faccia illuminata

function [MAT,ray_out] = reflection(MAT,ray,sat0)

                R.punti_imp_rif   = [];           
                R.face_inc_rif    = [];       % indice faccia illuminata da raggio
                R.ray_rif         = [];               % raggio riflesso


P_rif            = ray.p_source;    % punti partenza raggio (indice riga MAT.Centres)
rif              = ray.dir;         % direzione raggi

facce_source_idx = ray.face_source; % indice facce sorgenti  (sat0.Faces(idx))


for ii = 1:length(facce_source_idx)  % cicla per ogni raggio riflesso (delle facce ill)
    ss = 1;

    % indice riga delle facce la cui normale è < 0 con lanormale alla
    % faccia mesh e faccia posta davanti
    idx_ill = indici_facce_ill(rif(ii,:),sat0,P_rif(ii,:),facce_source_idx(ii));


    if ~isempty(idx_ill)

        % esclude le facce se non normale>0 e posizionate di fronte e le
        % ordina per distanza
        VV_rest_sort = Calcola_distanze(sat0,idx_ill,facce_source_idx(ii));
        idx_facce_pot = table2array(VV_rest_sort(:,2));

    end

       

    % verifica se il raggio colpisca la facce potenziali (ordinate da
    % quella più vicina)
    trovato = false; 

    if  exist('idx_facce_pot','var')
        if length(idx_facce_pot) > 1

            while trovato == false

                idx_riga_faccia     = VV_rest_sort.faccia_idx(ss);  % indice  riga della Faces illuminata dal rif
                idx_righe_Vertex    = sat0.Faces(idx_riga_faccia,:); % valori corrispondenti (indice riga Vertex)

                % coordinate mesh potenzialmente illuminata
                A        = sat0.Vertex(idx_righe_Vertex(1),:);
                B        = sat0.Vertex(idx_righe_Vertex(2),:);
                C        = sat0.Vertex(idx_righe_Vertex(3),:);

                % punto di intersezione raggio con piano mesh illumnato
                % potenziale

                if dot(rif(ii,:),sat0.Normals_mesh(idx_riga_faccia,:)) == 0
                    P_imp(1,:) = sat0.Centri_mesh(idx_riga_faccia,:);
                else
                    [P_imp(1,:),~] =  line_plane_intersection(rif(ii,:),P_rif(ii,:),...
                        sat0.Normals_mesh(idx_riga_faccia,:),A);
                end

                % verifica se il punto di impatto sta dentro la mesh
                % potenzialmnte illuminata
                if norm(P_imp) < 1e5  % escludi punti all'infinito
                    trovato = ver_inner(A,B,C,P_imp);
                end

                ss = ss+1;

                if ss > length(idx_facce_pot)
                    trovato = false;
                    break
                end

            end
        end

        clear idx_facce_pot

        if trovato == true

            % normale della faccia illuminata
            N        =  sat0.Normals_mesh(idx_riga_faccia,:);
            ray_rif  = 2*(dot(N,-rif(ii,:))).*N + rif(ii,:);

            if  dot(N,-rif(ii,:)) ~= 0  % Normale e raggio non complanari

                R.punti_imp_rif   = [R.punti_imp_rif ;P_imp];
                R.face_inc_rif    = [R.face_inc_rif  ;idx_riga_faccia];    % indice faccia illuminata da raggio
                R.ray_rif         = [R.ray_rif ;ray_rif];                  % raggio riflesso

            end

        end
    end
end

  
ray_out                 = c_Rays();
ray_out.face_source     = R.face_inc_rif    ;     % indice riga delle Face
ray_out.p_source        = R.punti_imp_rif ;       
ray_out.dir             = R.ray_rif  ;

MAT.face_source        = [MAT.source;R.face_inc_rif]; 
MAT.n_riflesso         = [MAT.n_riflesso; (MAT(end).n_riflesso(end)+1)*ones(length(R.face_inc_rif),1)]; 
end

