%Main Ray
clear;clc;
addpath('C:\Users\carlo\OneDrive - INAF - Istituto Nazionale di Astrofisica\G4S\Matlab_ray\Funzioni\')
addpath('C:\Users\carlo\OneDrive - INAF - Istituto Nazionale di Astrofisica\G4S\Matlab_ray\Dati Caomsol\')
% addpath('C:\Users\carlo\OneDrive - INAF - Istituto Nazionale di Astrofisica\G4S\Matlab_ray\Geometria\');

%% import material file
[MAT.Vertex ,MAT.Faces,MAT.matidx,Mat_tot]  = Read_mesh_MAT_txt('material_data_section.txt');

%% Matrice Sorgente raggi

ray_dir = [0 0 1];  % Line direction X Y(ala) Z
ray_start_x = linspace(0.4,0.4,10)';             % Point of the line
ray_start_y = linspace(-1.5,1.5,30)';            % Point of the line  (Wings)
ray_start_z = -1*ones(size(ray_start_x,2),1);    % Point of the line

% piano sprgente
% for tt = 1 : size(V_stl,1)
%  [ray_start_p(tt,:),~] = line_plane_intersection(ray_dir, V_stl(tt,:),[0 0 -1], [0 0 -1]); % punti proiezione dei Vertici su un piano
% end
% 
% face = [];
% tic
% RX = ray_start_p(:,1);
% for  p = 1:numel(RX)
%     ray_start = [ray_start_p(p,1),ray_start_p(p,2),ray_start_p(p,3)];
%     face_n  = Face_illuminated(F_stl,V_stl,ray_dir,ray_start);
%     if size(face_n,1)<1
%         break
%     else
%         face      = [face face_n];
%     end
% end


% % calcola faccia impatto e posizione
tic
face = [];
parfor  rx = 1:numel(ray_start_x)
      for ry = 1:numel(ray_start_y )
         for  rz = 1:numel(ray_start_z)
            ray_start = [ray_start_x(rx),ray_start_y(ry),ray_start_z(rz)];
            face_n  = Face_illuminated(Mat_tot,ray_dir,ray_start);

            if size(face_n,1)<1
                break
            else
             face      = [face face_n];
            end
        end
    end
end

toc

%% associa Materiale
% for jj = 1:length([face.area])
%     for kk = 1:size(Material,1)
%         dA = [face(jj).impact_P ];
%         dB = Material(kk,1:3);
%         L(kk) = norm(dA-dB);   % distanze del punto su face dal punto su Mater
%         [a, b]= (min(L));
%         face(jj).Mat = Material(b,4);
%     end
% end


% for jj = 1:length(V_stl)
%     for kk = 1:size(Material,1)
%         dA = V_stl(jj,1:3);
%         dB = Material(kk,1:3);
%         L(kk) = norm(dA-dB);   % distanze del punto su V_Stl dal punto su Mater
%         [a, b]= (min(L));
%         Mat(jj) = Material(b,4); % vettore Materiali con indici di V_stl
%     end
% end


% elimina facce interne
clear idx;

for kk = 1:length([face.area]) -1
    for jj = kk+1:length([face.area])
        if  sum(([face(kk).ray_source_P]) == ([face(jj).ray_source_P])) == 3
            if [face(kk).dist] < [face(jj).dist]
                idx(kk)  = jj;
            else
                idx(kk)  = kk;
            end
        end
    end
end

face_i = face;
if exist('idx','var')
    idx_u = unique(idx);
    face_i(idx_u(2:end)) = [];
end

%%

for fc = 1:size([face_i(1,:)],2)
    plot_tri([face_i(1,fc).vertex_A],[face_i(1,fc).vertex_B],[face_i(1,fc).vertex_C],'r',face_i(fc).view_fact)
end
axis equal

view(-ray_dir(1)/ray_dir(2) ,-rad2deg(asin(ray_dir(3)/v_leng(ray_dir))))


% figure(3)
% hold on
% for hh =1:length([face.area])
%     plot_tri_mat(face(hh))
% end

axis equal

