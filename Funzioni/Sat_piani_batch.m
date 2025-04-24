function Dist = Sat_piani_batch(A_ref, B_ref, C_ref, A_mat, B_mat, C_mat, P_mat)
    % Calcola i vettori normali per tutte le facce di riferimento
    v0_ref = B_ref - A_ref;
    v1_ref = C_ref - A_ref;
    normal_ref = cross(v0_ref, v1_ref, 2); % 1 x 3 vettore normale

    % Calcola i vettori normali per tutte le altre facce
    v0 = B_mat - A_mat;
    v1 = C_mat - A_mat;
    normal = cross(v0, v1, 2); % m x 3 vettori normali

    % Normalizza i vettori normali
    normal_ref = normal_ref / norm(normal_ref);
    normal_norm = sqrt(sum(normal.^2, 2));
    normal = normal ./ normal_norm;

    % Calcola le distanze usando il prodotto scalare tra i normali e i punti
    Dist = dot(normal, repmat(normal_ref, size(normal, 1), 1), 2);
end
