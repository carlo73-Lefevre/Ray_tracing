function IDX_Facce_ill_vect = dentro_list_vect_002(A, B, C, P,Sat)

%% Passo 1: Calcola i vettori dei triangoli
v0 = C - A; % n x 3  
v1 = B - A; % n x 3  

% Passo 2: Prepara i punti P
m = size(P, 1); 
n = size(A, 1);

% Passo 3: Calcola v2 usando broadcasting esplicito
% A_expanded avrà dimensione n x 1 x 3, mentre P avrà dimensione 1 x m x 3
A_expanded = reshape(A, [n, 1, 3]);
P_expanded = reshape(P, [1, m, 3]);

% Sottrazione vettoriale con broadcasting per ottenere v2 (n x m x 3)
v2 = P_expanded - A_expanded;

% Passo 4: Calcola i prodotti scalari
% Espandi v0 e v1 per farli combaciare con v2
v0_expanded = reshape(v0, [n, 1, 3]);  % v0 diventa n x 1 x 3
v1_expanded = reshape(v1, [n, 1, 3]);  % v1 diventa n x 1 x 3

% Prodotti scalari tra v0 e v2, v1 e v2 (dot02, dot12 avranno dimensione n x m)
dot00 = sum(v0_expanded .* v0_expanded, 3); % n x 1  
dot01 = sum(v0_expanded .* v1_expanded, 3);  
dot11 = sum(v1_expanded .* v1_expanded, 3);  
dot02 = sum(v0_expanded .* v2, 3); % n x m  
dot12 = sum(v1_expanded .* v2, 3);

% Passo 5: Calcola le coordinate baricentriche
denom = dot00 .* dot11 - dot01 .* dot01;
u = (dot11 .* dot02 - dot01 .* dot12) ./ denom;
v = (dot00 .* dot12 - dot01 .* dot02) ./ denom;

% Leva le sovrapposizioni
u(1:n+1:end) = -1;
v(1:n+1:end) = -1;

% Passo 6: Determina se i punti sono all'interno dei triangoli
inside = (u > 0) & (v > 0) & (u + v < 1);

% Passo 7: Identifica i punti all'interno di almeno un triangolo
points_inside_any_triangle = any(inside, 1);
indices_of_points_inside = find(points_inside_any_triangle);

% Associazioni tra triangoli e punti
[triangle_indices, point_indices] = find(inside);

% Verifica che il triangolo in esame sia davanti alla faccia oscurata
idx_pt = triangle_indices < point_indices;
Tri_ok = triangle_indices;
Pt_ok  = point_indices;

Tri_ok(idx_pt) = Pt_ok(idx_pt);

IDX_Facce_ill_vect = Tri_ok;

end
