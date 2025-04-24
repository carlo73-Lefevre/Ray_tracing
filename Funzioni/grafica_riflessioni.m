function p1 = grafica_riflessioni(Sat,ray1)



p0 = patch('Faces',Sat.Faces,'Vertices',Sat.Vertex); %,'EdgeColor','none'
p0.FaceColor = [0.7 0.7 0.7]; hold on; grid on;axis equal;
p1 = patch('Faces',Sat.Faces(ray1.face_source,:),'Vertices',Sat.Vertex);
p1.FaceColor = 'y'; hold on; grid on;axis equal;

% plot_quiv(ray1.p_source,ray1.dir,1,'r');


end