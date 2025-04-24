function IDX_Facce_ill_vect = dentro_list_vect(Sat0,IDX_Facce_ill,A,B,C,P)
% Supponiamo che A, B, C, e P siano definiti come descritto  

%% Passo 1: Calcola i vettori dei triangoli  

v0 = C - A; % n x 3  
v1 = B - A; % n x 3  
% v2 = A - B;

% Passo 2: Prepara i punti P  
% P = A; % m x 3  

m = size(P, 1);  
n = size(A,1);

% Passo 3: Calcola v2  
A_expanded = reshape(A, [n, 1, 3]);  
A_expanded = repmat(A_expanded, [1, m, 1]);  

P_expanded = reshape(P, [1, m, 3]);  
P_expanded = repmat(P_expanded, [n, 1, 1]);  

v2 = P_expanded - A_expanded; % n x m x 3  

% Passo 4: Calcola i prodotti scalari  
v0_expanded = repmat(reshape(v0, [n, 1, 3]), [1, m, 1]);  
v1_expanded = repmat(reshape(v1, [n, 1, 3]), [1, m, 1]);  

dot00 = sum(v0_expanded .* v0_expanded, 3); % n x m  
dot01 = sum(v0_expanded .* v1_expanded, 3);  
dot11 = sum(v1_expanded .* v1_expanded, 3);  
dot02 = sum(v0_expanded .* v2, 3);  
dot12 = sum(v1_expanded .* v2, 3);  

% Passo 5: Calcola le coordinate baricentriche do P 
denom = dot00 .* dot11 - dot01 .* dot01;  
u = (dot11 .* dot02 - dot01 .* dot12) ./ denom;  
v = (dot00 .* dot12 - dot01 .* dot02) ./ denom;  

% leva le sovrapposizioni
mu = size(u,1);
u(1:mu+1:end) = -1;
v(1:mu+1:end) = -1;


% Passo 6: Determina se i punti sono all'interno dei triangoli  
inside = (u > 0 ) & (v > 0) & (u + v < 1); % n x m  

% Passo 7: Identifica i punti all'interno di almeno un triangolo  
points_inside_any_triangle = any(inside, 1); % 1 x m  
indices_of_points_inside = find(points_inside_any_triangle);  

% Associazioni tra triangoli e punti  
[triangle_indices, point_indices] = find(inside);  % triamgoli che mettono in ombra i punti
% idx_ext = (triangle_indices- point_indices~=0);



% verifica che il triangolo in esame sia davanti alla faccia osurata, in
% caso contrario scambiali
idx_pt = triangle_indices < point_indices;
Tri_ok = triangle_indices;
Pt_ok  = point_indices;

Tri_ok(idx_pt) = Pt_ok(idx_pt);





IDX_Facce_ill_vect = Tri_ok;



%%
end
