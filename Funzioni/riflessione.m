function riflessi_pot = riflessione(MAT_old)

% verifica se il raggio riflesso impatta
for ii = 1 : size(MAT_old.facce_restanti,2)
    direzione_rif = MAT_old.ray_rif(ii,:);
    piano_proiezione = direzione_rif*10;
    MAT_spare = facce_illuminate_potenziali(MAT_old,direzione_rif);

end

riflessi_pot = MAT_spare.facce_potenzialmente_illuminate;

end