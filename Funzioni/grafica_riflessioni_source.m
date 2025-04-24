function p1 = grafica_riflessioni_source(Sat,ray1,ray2)



p0 = patch('Faces',Sat.Faces,'Vertices',Sat.Vertex); %,'EdgeColor','none'
p0.FaceColor = [0.7 0.7 0.7]; hold on; grid on;axis equal;

p1 = patch('Faces',Sat.Faces(ray2.face_source,:),'Vertices',Sat.Vertex); %faccia di arrivo
p1.FaceColor = 'r'; hold on; grid on;axis equal;
%% calcola la direzione del raggio di provenienza
reflect_vector = reflect(Sat.Normals_mesh(ray2.face_source,:),ray2.dir);

%%
plot_quiv(Sat.Centri_mesh(ray2.face_from,:),reflect_vector,3);



function reflect_vector = reflect(surface_normal, incident_vector)

    reflect_vector = zeros(size(surface_normal,1),3);

for j = 1:size(surface_normal,1)
    reflect_vector(j,:) = incident_vector(j,:) - 2 .* (incident_vector(j,:) * surface_normal(j,:)') .* surface_normal(j,:);
end

end


end