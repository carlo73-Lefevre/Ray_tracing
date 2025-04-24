function intersection = Intersection_CGPT(u, N, n, M)
    % Calcola il parametro t per l'intersezione
    t = dot(n, (M - N)) / dot(n, u);

    % Calcola le coordinate del punto di intersezione
    intersection = N + t * u;
end