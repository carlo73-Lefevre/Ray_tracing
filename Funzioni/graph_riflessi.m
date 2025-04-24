
function graph_riflessi(Sat,MAT,ray)

figure
clf
hold on;
 p = patch('Faces',Sat.Faces,'Vertices',Sat.Vertex,'EdgeColor','none');%
 p.FaceColor =  'k';
%  p.FaceAlpha = 0.8;
% facce illuminate dai riflessi
   % facce illuminate
    p5 = patch('Faces',Sat.Faces(MAT.facce_restanti,:),'Vertices',Sat.Vertex); %,'EdgeColor','none'
    p5.FaceColor = [0.7,0.7,0.7];
    p5.FaceAlpha = 0.9;

    plot_quiv(ray.p_source,ray.dir,1)  % raggi riflessi 1
%      plot_tri(rif3.p_imp)

box on;
axis equal
legend('facce in ombra','facce iluminate dirette','facce raggi riflessi 1')
view(-120,28)
end