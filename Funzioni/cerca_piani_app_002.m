function App = cerca_piani_app_002(idx_faces, center_n, App, f, Sat)

% Preleva i vertici e i centri del satellite
A = Sat.Vertex(Sat.Faces(:,1),:);
B = Sat.Vertex(Sat.Faces(:,2),:);
C = Sat.Vertex(Sat.Faces(:,3),:);
P = Sat.Centers_mesh;

% Condizione di arresto: se f è maggiore o uguale alla lunghezza di center_n o se tutte le facce sono state marcate
if f >= length(center_n) || all(idx_faces)
    return; % Fermati se tutte le facce sono state processate o f ha superato la lunghezza
end

% Seleziona solo le facce non ancora marcate
Fac = find(idx_faces == 0);

if ~isempty(Fac)
    % Prendi le coordinate del primo triangolo di riferimento
    A_ref = A(Fac(1), :);
    B_ref = B(Fac(1), :);
    C_ref = C(Fac(1), :);

    % Prepara le matrici per il calcolo vettoriale
    A_mat = A(Fac, :); % m x 3
    B_mat = B(Fac, :); % m x 3
    C_mat = C(Fac, :); % m x 3
    P_mat = P(Fac, :); % m x 3

    % Calcola le distanze per tutte le facce in modo vettoriale
    Dist = Sat_piani_batch(A_ref, B_ref, C_ref, A_mat, B_mat, C_mat, P_mat);

    % Trova le facce complanari con tolleranza
    idx_c = abs(Dist) < 1e-14;

    % Assegna il primo indice delle facce complanari ad `App`
    if any(idx_c)
        App(Fac(idx_c)) = Fac(1); % Usa l'indice del primo triangolo non marcato
    end

    % Aggiorna `idx_faces` per marcare le facce complanari
    idx_faces(Fac(idx_c)) = true;

    % Richiama ricorsivamente la funzione con il nuovo indice
    % Solo se ci sono ancora facce non marcate
    % if any(~idx_faces)
    %     App = cerca_piani_app_002(idx_faces, center_n, App, sum(idx_faces), Sat);
    % end
end

end