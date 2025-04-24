%% Calcola le facce che non sono in ombra.
%  La funzione aggiorna l'idice delle facce illuminati. Queste sono quelle
%  che restano dal totale

function idx_rest = escludi_facce(MAT,nray,cont)
% cont : contatore delle facce analizzate. Se alla chiamata precedente non
% è stata trovata nessuna faccia in ombra, il cont=cont+1 

IDX_Facce_ill = MAT.facce_restanti;       % variabile dinamica, indici facce restanti

% proiezioni del triancolino fem sul piano sorgente
Faccia_idx   = IDX_Facce_ill(cont);   % questo indice è lo steso pe le Facce per i ray_dir_v e ray_start_p
raggio_dir   = nray.dir_rif1(Faccia_idx,:);             %  MAT.ray_dir_v(Faccia_idx,:);
raggio_start = nray.p_source_rif1(Faccia_idx,:);        %  MAT.ray_start_p(Faccia_idx,:);
idx_P        = MAT.Faces(Faccia_idx,:);

% posizione Centro faccia sorgente
C0   = MAT.Centri_mesh(Faccia_idx,:);


for jj = 1:length(MAT.Faces)  % sulle facce illuminate potenzialmente verifica se il triangolo contiene punti sorgente
    idx_F = MAT.Faces(jj,:);
    % vertici delle facce
    Va = MAT.Vertex(idx_F(1),:);
    Vb = MAT.Vertex(idx_F(2),:);
    Vc = MAT.Vertex(idx_F(3),:);

    Dentro(jj) = ver_inner(Va,Vb,Vc,C0); % verifica se la facia è colpita dal rggio riflesso

    if MAT.Faces(jj) == Faccia_idx
        Dentro(jj) = false;
    end

end

[Den_val,Den_pos]= find(double(Dentro) == 0); % indici di idx_fc_ill_n da eliminare (Den_b)

if Den_val == true
    Dentro_1 =  Den_pos;
end

if exist("Dentro_1",'var')
 % cancella tutte le facce in ombra e rimani in cont attuale
    MAT.facce_restanti(Dentro_1) = [];

else
 % analizza la faccia illuminata successiva delle rimanenti visto che quella attuale non è da elimnare
    cont = cont+1;

end

idx_rest = MAT.facce_restanti;

if cont < length(MAT.facce_restanti)
    idx_rest = escludi_facce(MAT,nray,cont);
end

end