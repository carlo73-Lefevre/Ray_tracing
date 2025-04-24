function ray_rif = trova_riflessi(sat,ray_dir)

for Fnum = 1:size(sat.Faces,1) % cicla per tutte le facce

    N =  sat.Normals_mesh(Fnum,:);

    ray_rif(Fnum,:)  = 2*(dot(N,-ray_dir(Fnum,:))).*N + ray_dir(Fnum,:); % formul della riflessione (con direzione uscente)
    ray_rif(Fnum,:)  = ray_rif(Fnum,:) / norm(ray_rif(Fnum,:));  % angolo raggio riflesso

end
end