function ray_rif = raggi_riflessi(CNorm,ray_dir)
% Calcola le direzioni di riflessione
N = CNorm;
ray_rif = 2 * dot(N, -ray_dir, 2) .* N + ray_dir;
ray_rif = ray_rif ./ vecnorm(ray_rif, 2, 2);

end