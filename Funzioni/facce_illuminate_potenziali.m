function [MAT,loc] = facce_illuminate_potenziali(sat,ray)
%%%
            
    DDot = sum(sat.Normals_mesh'.*( ray.dir'./norm(ray.dir))); % Prodotto scalare facce sorgenti
    idx_fc_ill = find(DDot < -1e-5); % Indici facce potenzialmente illuminate

    % Calcola le distanze tra i centri delle facce e la sorgente di luce
    D = vecnorm(ray.p_source - sat.Centers_mesh, 2, 2);

    % Ordina le distanze e ottieni gli indici
    [sorted_D, sorted_idx] = sort(D);

    % Trova le facce potenziali tra quelle rivolte ne verso giusto
    [~, loc] = ismember(sorted_idx, idx_fc_ill);

    % Filtra solo le facce potenzialmente illuminate
    facce_ill_idx  = sorted_idx(loc > 0);
    facce_restanti = sorted_idx(loc == 0);

    % Crea la strutt ura di output
    MAT.Dis_sort_ill        = facce_ill_idx;
    MAT.facce_restanti      = facce_restanti;
    MAT.ordine_per_distanza = sorted_idx(loc > 0);



end