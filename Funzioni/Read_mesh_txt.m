function [Vertex,Faces] = Read_mesh_txt(filename)
%% test stl import
T = readtable(filename);
L = table2array(T(2,3));

mesh = load(filename);%test_mesh.txt
Vertex = mesh(1:L ,1:3);
Faces =  mesh(L +1:end,1:3);
% triangoli = [];
% for k=1:size(Faces,1)
%     for j=1:3
%         triangoli = [triangoli; Vertex(Faces(k,j),1),Vertex(Faces(k,j),2),Vertex(Faces(k,j),3)];
%     end
% end

fig = figure(1);
clf;
ax = axes('Parent',fig);
% fill3(triangoli(:,1),triangoli(:,2),triangoli(:,3),'r')
p = patch('Faces', Faces, 'Vertices', Vertex);
p.Parent = ax;
hold(ax,'on');
p.FaceColor = [0.9 0.95 0.95];
p.FaceAlpha = 0.8;

view([-20 -45])
axis equal
lighting gouraud
title('Mesh')
grid on