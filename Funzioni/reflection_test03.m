
%% Calcola raggio riflesso e indice faccia illuminata

function ray_out = reflection_test03 (ray,sat0,varargin)

R.p_source      = [];
R.face_source   = [];   % indice faccia illuminata da raggio
R.dir           = [];   % raggio riflesso
R.face_from     = [];   % faccia di provenienza


P_rif            = ray.p_source;    % punti partenza raggio (indice riga MAT.Centres)
rif_dir          = ray.dir;         % direzione raggi

facce_source_idx = ray.face_source; % indice facce sorgenti  (sat0.Faces(idx))

for ii = 1:length(facce_source_idx)  % cicla per ogni raggio riflesso (delle facce ill)


    % indice riga delle facce la cui normale è < 0 con la normale alla
    % faccia mesh e faccia posta davanti

    idx_ill = indici_facce_ill(rif_dir(ii,:),sat0,P_rif(ii,:),facce_source_idx(ii));


    %% proietta tutti i centri delle mesh delle facce idx_ill su di un piano
    % lontano e di fronte a rif(ii,:)
    if ~isempty(idx_ill)
        % esclude le facce se non normale>0 e posizionate di fronte e le
        % ordina per distanza
        VV_rest_sort  = Calcola_distanze(sat0,idx_ill,facce_source_idx(ii));
        % idx_facce_pot = table2array(VV_rest_sort(:,2));
                idx_facce_pot = (VV_rest_sort(:,2));
    else 
        continue
    end



    % verifica se il raggio colpisca la facce potenziali (ordinate da
    % quella più vicina)
    trovato = false;

    if  exist('idx_facce_pot','var')

        if length(idx_facce_pot) > 1
            ss = 1;  %contatore delle facce potenziali analizzate

            while trovato == false

                idx_riga_faccia     = idx_facce_pot(ss);  % indice  riga della Faces illuminata dal rif
                idx_righe_Vertex    = sat0.Faces(idx_riga_faccia,:); % valori corrispondenti (indice riga Vertex)

                % coordinate mesh potenzialmente illuminata 
                A        = sat0.Vertex(idx_righe_Vertex(1),:);
                B        = sat0.Vertex(idx_righe_Vertex(2),:);
                C        = sat0.Vertex(idx_righe_Vertex(3),:);

                % punto di intersezione raggio con piano mesh illumnato
                % potenziale
                   
                % se sono ortigonale raggio e normale 
                if abs(dot(rif_dir(ii,:),sat0.Normals_mesh(idx_riga_faccia,:))) < 1e-8 %dot(rif(ii,:),sat0.Normals_mesh(idx_riga_faccia,:)) == 0
                    
                      trovato = false;
                      ss = ss+1;
                      break
%                     P_imp(1,:) = sat0.Centri_mesh(idx_riga_faccia,:);
                else  % se non sono ortogonale trova il unto di impatto
                    % [P_imp(1,:),~] =  line_plane_intersection(rif(ii,:),P_rif(ii,:),...
                    %     sat0.Normals_mesh(idx_riga_faccia,:),A);
                    P_imp(1,:) =   Intersection_CGPT(rif_dir(ii,:),P_rif(ii,:),...
                                   sat0.Normals_mesh(idx_riga_faccia,:),A);
                end

                % verifica se il punto di impatto sta dentro la mesh
                % potenzialmnte illuminata
                if norm(P_imp) < 1e5  % escludi punti all'infinito
                    trovato = in_or_out(A,B,C,P_imp);
                end

                ss = ss+1;

                if ss > length(idx_facce_pot) || dot(rif_dir(ii,:),sat0.Normals_mesh(idx_riga_faccia,:)) == 0
                    trovato = false;
                    break
                end

            end
        end

        clear idx_facce_pot
        % una volta trovato quale facci acolpishe il raggio ii, si deve
        % verificare quali centri faccia rientrano nella proiezione della
        % mesh del raggio sorgente sul piano complanare alla faccia colpita

        if trovato == true  % se P_imp è dentro la mesh analizzata

            % Proietta tutti i centri delle facce complanari alla faccia
            % illuminata sul piano del raggio sorgente

            % indici delle facce complanari alla colpita
            Compl_facce_idx = sat0.Face_coplanar == sat0.Face_coplanar(idx_riga_faccia);
            idx_compl = find( Compl_facce_idx == 1);  % indice facce complanari alla colpita

            Rest.facce_restanti = idx_compl;  % facce complanari da analizzzare

            N_ray_mesh                =   sat0.Normals_mesh(ray.face_source(ii),:);  % normal alla facccia del raggio
            Centri_mesh               =   ones(size(sat0.Faces,1),3) .* [0 0 1000];
            Centri_mesh(idx_compl,:)  =   sat0.Centers_mesh(idx_compl,:); % Centri delle complanari
            P_source                  =   P_rif(ii,:);
            ray_dir                   =   rif_dir(ii,:);
            ray_face                  =   ray.face_source(ii);

            % verifica che la faccia non sia stata già illuminata da
            % un'altro raggio
            if size(R.face_source,2)>0
                idx_out = [];
                for aa = 1 : size(R.face_source,1) % cicla per le facce complanari
                    idx = find(Rest.facce_restanti == R.face_source(aa));
                    if size(idx,1)>0

                        idx_out(end+1) =  find(Rest.facce_restanti == R.face_source(aa));% indici già analizzati
                    
                    end
                end

                Rest.facce_restanti(idx_out) = [];

            end


            if size(Rest.facce_restanti,1)>0 && size(Rest.facce_restanti,2)>0
                %%%%%%% veriifica quali facce (indice) del piano complanare sono
                % nella mesh sorgente proiettata %%%%%%%%

                idx_face_rif = escludi_facce_rif(sat0,Rest,Centri_mesh,N_ray_mesh,P_source,ray_dir,ray_face,varargin);
   
   
                % se non ne trova salva solo la prima normale della faccia illuminata
                if isempty(idx_face_rif)
                    idx_face_rif = idx_riga_faccia;
                end
      
                N            =  sat0.Normals_mesh(idx_face_rif(1),:);  % serve solo la prima perchè sono tutte uguali
                ray_rif_dir  = 2*(dot(N,-rif_dir(ii,:))).*N + rif_dir(ii,:);
                ray_rif      = ones(length(idx_face_rif),3).* ray_rif_dir;

               %  centro della mesh colpita (delle complanari analizzate)
                P_imp =  sat0.Centers_mesh(idx_face_rif,:);

                % registra le facce di provenienza
                idx_facce_provenienaa = facce_source_idx(ii) * ones(size(P_imp,1),1); 

                if  dot(N,-rif_dir(ii,:)) ~= 0  % Normale e raggio non complanari

                    R.p_source      = [R.p_source ;P_imp];
                    R.face_source   = [R.face_source  ;idx_face_rif'];        % indice faccia illuminata da raggio
                    R.dir           = [R.dir ;ray_rif];                       % raggio riflesso
                    R.face_from     = [R.face_from  ;idx_facce_provenienaa];  % faccia di provenienza

                    clear idx_face_rif facce_source ray_rif P_imp

                end
            end
        end
    end
end


ray_out                 = c_Rays();
ray_out.face_source     = R.face_source    ;     % indice riga delle Face
ray_out.p_source        = R.p_source ;
ray_out.dir             = R.dir  ;
ray_out.face_from       = R.face_from  ;

end

