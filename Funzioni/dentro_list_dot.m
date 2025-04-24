function Sat0_ill = dentro_list_dot(Sat0, n_dir, rays_source, cont, varargin)
% P0: posizione di partenza del raggio sorgente sul piano sorgente
% IDX_Facce_ill: indici delle facce di Sat0 illuminate

% Recupero delle facce illuminate e i loro vertici
% Facce_ill_vert = Sat0.Faces(IDX_Facce_ill,:);

% Estrazione dei punti della mesh relativi ai vertici A, B, C
% A = Sat0.Vertex(Facce_ill_vert(:,1),:);
% B = Sat0.Vertex(Facce_ill_vert(:,2),:);
% C = Sat0.Vertex(Facce_ill_vert(:,3),:);

A = Sat0.Vertex(Sat0.Faces(:,1),:);
B = Sat0.Vertex(Sat0.Faces(:,2),:);
C = Sat0.Vertex(Sat0.Faces(:,3),:);


cont
% Calcolo della distanza lungo i vettori dei raggi (proiezione dei vertici)
t_F_A = sum(n_dir .* (rays_source - A), 2) ./ sum(n_dir .* n_dir, 2);
t_F_B = sum(n_dir .* (rays_source - B), 2) ./ sum(n_dir .* n_dir, 2);
t_F_C = sum(n_dir .* (rays_source - C), 2) ./ sum(n_dir .* n_dir, 2);

% Calcolo dei nuovi vertici sul piano sorgente
A = A + t_F_A .* n_dir;
B = B + t_F_B .* n_dir;
C = C + t_F_C .* n_dir;

% Inizializzazione variabili di output
% IDX_Facce_ill_n = IDX_Facce_ill;  % Copia iniziale degli indici
 % Sat0.Faces_n =  Sat0.Faces;
c = cont; 
% Ciclo per verificare le facce in ombra prodotto dalle mesh illuminate
while c > 1

    At = A(c-1,:);
    Bt = B(c-1,:);
    Ct = C(c-1,:);

    % vettori edge triangoli
    ab = Bt - At;
    bc = Ct - Bt;
    ca = At - Ct;


    % Calcolo dei vettori tra il triangolo in esame (c) e tutti gli altri triangoli

    P1 = A(1:c,:);
    P2 = B(1:c,:);
    P3 = C(1:c,:);

    % primo punto del triangolo
    vaP1 =  At - P1 ;
    vbP1 =  Bt - P1 ;
    vcP1 =  Ct - P1 ;

    abr =  repmat(ab,size(vaP1,1), 1);
    bcr =  repmat(bc,size(vbP1,1), 1);
    car =  repmat(ca,size(vcP1,1), 1);

    ScalPa1 = sum(n_dir(1:size(vaP1,1),:) .* cross(abr,vaP1), 2);
    ScalPb1 = sum(n_dir(1:size(vaP1,1),:) .* cross(bcr,vbP1), 2);
    ScalPc1 = sum(n_dir(1:size(vaP1,1),:) .* cross(car,vcP1), 2);

    % secondo punto del triangolo
    vaP2 =  At - P2 ;
    vbP2 =  Bt - P2 ;
    vcP2 =  Ct - P2 ;

    ScalPa2 = sum(n_dir(1:size(vaP1,1),:)  .* cross(abr,vaP2), 2);
    ScalPb2 = sum(n_dir(1:size(vaP1,1),:)  .* cross(bcr,vbP2), 2);
    ScalPc2 = sum(n_dir(1:size(vaP1,1),:)  .* cross(car,vcP2), 2);

    % terzo punto del triangolo
    vaP3 =  At - P3;
    vbP3 =  Bt - P3;
    vcP3 =  Ct - P3;

    ScalPa3 = sum(n_dir(1:size(vaP1,1),:)  .* cross(abr,vaP3), 2);
    ScalPb3 = sum(n_dir(1:size(vaP1,1),:)  .* cross(bcr,vbP3), 2);
    ScalPc3 = sum(n_dir(1:size(vaP1,1),:)  .* cross(car,vcP3), 2);


    % Verifica se i vertici sono interni (tutti i prodotti scalari positivi)
    idx_a = all([ScalPa1 ScalPb1 ScalPc1] >  1e-4, 2);
    idx_b = all([ScalPa2 ScalPb2 ScalPc2] >  1e-4, 2);
    idx_c = all([ScalPa3 ScalPb3 ScalPc3] >  1e-4, 2);
    % Calcolo dell'ombreggiatura combinando i risultati
    idx_shadow = or(or(idx_a, idx_b), idx_c);  % indici della matrice Sat.Faces in ombra
    % Aggiorna gli indici delle facce illuminate (elimina le facce in ombra)

    % idx_F = Sat0.Vertex(idx_shadow,:); % indice riga della Delle facce in ombra (in Sat0.Face)

    if sum(idx_shadow>0)
        % [idf1,idf11] = find(Sat0.Faces(:,1) == idx_F(1,:) );
        % [idf2,idf22] = find(Sat0.Faces(:,2) == idx_F(2,:) );
        % [idf3,idf33] = find(Sat0.Faces(:,3) == idx_F(3,:) );

        % idx_IDX = [idf1;idf2;idf3];

   
        % rays_source( IDX_Facce_ill_n(unique(idx_IDX)),:)      = [];
        
        n_dir(1:sum(idx_shadow),:)     = [];
        rays_source(idx_shadow,:)      = [];
        % IDX_Facce_ill(idx_shadow)     = [];

        Sat0.Faces(idx_shadow,:)      = [];
        

        c = size(Sat0.Faces,1);
        Sat0 = dentro_list_dot(Sat0, n_dir, rays_source, c, varargin)

    else
        c = c-1;        
        Sat0 = dentro_list_dot(Sat0, n_dir, rays_source, c, varargin)

    end



    % end
    % jj = jj+1;
end

%Visualizzazione delle facce




% % Aggiorna A, B, C per la nuova lista di facce illuminate
% A = A(~idx_shadow,:);
% B = B(~idx_shadow,:);
% C = C(~idx_shadow,:);
%
% % Ricalcola i lati del triangolo per la nuova lista di facce illuminate
% ab = B - A;
% bc = C - B;
% ca = A - C;


% IDX_Facce_ill_end = IDX_Facce_ill_n;

end
