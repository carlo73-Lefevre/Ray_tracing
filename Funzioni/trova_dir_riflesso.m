function [CNorm,ray_rif] = trova_dir_riflesso(sat,ray_dir)

% Estrai i nodi delle facce
A = sat.Vertex(sat.Faces(:, 1), :);
B = sat.Vertex(sat.Faces(:, 2), :);
C = sat.Vertex(sat.Faces(:, 3), :);

% Calcola i vettori rappresentanti i lati delle facce
v1 = vector(A, B);
v2 = vector(B, C);
v3 = vector(C, A);

% Calcola le normali alle facce
CNorm = -cross(v1, v3) ./ (vecnorm(v1, 2, 2) .* vecnorm(v2, 2, 2));
CNorm = CNorm ./ vecnorm(CNorm, 2, 2);

% Calcola le direzioni di riflessione
N = CNorm;
ray_rif = 2 * dot(N, -ray_dir, 2) .* N + ray_dir;
ray_rif = ray_rif ./ vecnorm(ray_rif, 2, 2);


% if Fnum > 700
    % figure(1)
    % clf
    % hold on
    % p0 = patch('Faces',sat.Faces,'Vertices',sat.Vertex);
    % p0.FaceColor = [0.7 0.7 0.7];
    % p0.FaceAlpha = 0.8;
    % % p1 = patch('Faces',sat.Faces(Fnum,:),'Vertices',sat.Vertex);
    % % p1.FaceColor = 'r';
    % plot_quiv(sat.Centri_mesh,CNorm,2)
    % view(-120,28)
    % box on;
    % axis equal
% end

end

