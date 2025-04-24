%%
face =  7;

figure(1)
clf;
p0 = patch('Faces',Sat.Faces(4072:end,:),'Vertices',Sat.Vertex); % facce colpite dai riflessi
p0.FaceColor = [0.7,0.7,0.7];
p0.FaceAlpha = 1;
hold on; grid on;

p1 = patch('Faces',Sat.Faces(facce_source_idx(ii),:),'Vertices',Sat.Vertex); % facce colpite dai riflessi
p1.FaceColor = [0.2, 0.2, 1];
axis equal;
view(-60, 45);

plot_quiv(ray.p_source(ii,:),rif_dir(ii,:))

 plot3(P_rif(ii,1),P_rif(ii,2),P_rif(ii,3),'o')




% p3 = patch('Faces', Sat.Faces(idx_facce_pot(Dentro_idx_A(2)),:),'Vertices',Sat.Vertex); % facce colpite dai riflessi
% p3.FaceColor = 'r';
axis equal;
view(-60, 45);


plot_quiv(Sat.Centers_mesh(idx_facce_pot(Dentro_idx_A(1)),:),dir_ref)

p4 = patch('Faces', Sat.Faces(idx_facce_pot(Dentro_idx_A(1)),:),'Vertices',Sat.Vertex); % facce colpite dai riflessi
p4.FaceColor = 'r';

%%
figure(1)
clf
plot3(A(Tri_ok(1),1),A(Tri_ok(1),2),A(Tri_ok(1),3),'+')
hold on;
plot3(B(Tri_ok(1),1),B(Tri_ok(1),2),B(Tri_ok(1),3),'+')
plot3(C(Tri_ok(1),1),C(Tri_ok(1),2),C(Tri_ok(1),3),'+')
plot3(P(:,1),P(:,2),P(:,3),'o')

%% seconda riflessione
clf
cc=54;
p4 = patch('Faces', Sat_b(1,2).Faces,'Vertices',Sat.Vertex); % facce colpite dai riflessi
p4.FaceColor = [0.8 0.7 0.7];
hold on
p5 = patch('Faces', Sat_b(1,2).Faces(ray2(1,cc).face_from(:),:),'Vertices',Sat.Vertex); % facce colpite dai riflessi
p5.FaceColor = 'y';

view(-60, 45);
p5 = patch('Faces', Sat_b(1,2).Faces(ray2(1,1).face_source(cc),:),'Vertices',Sat.Vertex); % facce colpite dai riflessi
p5.FaceColor = 'r';
axis equal;


plot_quiv(ray2(1,cc).p_from,ray2(1,cc).dir_from)
% plot_quiv(ray2(1,1).p_source(cc,:),ray2(1,1).dir(cc,:),'r')