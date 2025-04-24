% calcola accelerazioni
function [ACC_S,ACC_N,Area_tot] = SRP(Sat,ray0,ray1,sun_dir,FI_mod,M)
% ray 0 raggi di partenza per calcolare il cos(beta)
% ray 1 raggi riflessi per avere le facce di from(da dove)

% const   = read_par('constant_SI.par');              % carica le costanti SI
% astrody = read_par('astrody_SI.par');    

FI = 1360.8;        % ;  % W/m^2
c_const = physconst('LightSpeed');  % [m/s]


% vertici facce di impatto raggio ray2 (ossia da dove partirà)
a = Sat.Vertex(Sat.Faces(ray1.face_source,1),:);
b = Sat.Vertex(Sat.Faces(ray1.face_source,2),:);
c = Sat.Vertex(Sat.Faces(ray1.face_source,3),:);

Acc_s = [];
Area_tot = [];

% figure(110)
% clf;
% p0 = patch('Faces',Sat.Faces,'Vertices',Sat.Vertex); %,'EdgeColor','none'
% p0.FaceColor = [0.7 0.7 0.7]; hold on; grid on;axis equal;
%  hold on;
% p1 = patch('Faces',Sat.Faces(ray1.face_source,:),'Vertices',Sat.Vertex); %,'EdgeColor','none'
% p1.FaceColor = [0.2 0.3 0.7];
% 
%  X = [1 0 0];Y =[0 1 0]; Z=[0 0 1];
% quiver3([0 0 0],[0 0 0],[0 0 0],X,Y,Z,6);
% title('facce illuminate in SRP funzione')

for jj =  1 : length(a)

    n       = Sat.Normals_mesh(ray1.face_source(jj),:); 
    Area    = calc_area(a(jj,:),b(jj,:),c(jj,:));
    mat_idx = Sat.materials_idx(Sat.Faces(ray1.face_source(jj),1));
    alpha   = Sat.optical_par(mat_idx).alfa;
    rho     = Sat.optical_par(mat_idx).rho;
    delta   = Sat.optical_par(mat_idx).delta;
    cos_b   = dot(n,-ray0.dir(jj,:) ); % angolo tra N e il vettore satellite_sole (o raggio incidente)
    

    Acc_n.S(jj)       = (-FI*FI_mod/(M*c_const)) * (alpha+delta) *Area * cos_b;                   % SRP sun 
    Acc_n.N(jj,:)     = (-FI*FI_mod/(M*c_const)) * 2*(delta/3 +rho*cos_b)*Area*abs(cos_b)*n;      % SRP N (delle Aree)  *FI_mod
    Acc_n.mat(jj)     = mat_idx;
    Acc_n.A(jj)       = Area; 
    Acc_n.cos_b(jj)   = cos_b;
    Area_tot          = sum([Area_tot Area]);



% p0 = patch('Faces',Sat.Faces,'Vertices',Sat.Vertex);
% p0.FaceColor = [0.5 0.5 0.5]; hold on; grid on;axis equal;
% 
% p1 = patch('Faces',Sat.Faces(ray1.face_source(ii),:),'Vertices',Sat.Vertex);
% p1.FaceColor = [mod(randn(1),1) mod(randn(1),1) mod(randn(1),1)]; 
% plot_quiv(ray1.p_source(ii,:),ray1.dir(ii,:),3);
% view(-61,12)


end

Acc_s_X = sum(Acc_n.S)*-sun_dir(1);
Acc_s_Y = sum(Acc_n.S)*-sun_dir(2);
Acc_s_Z = sum(Acc_n.S)*-sun_dir(3);

ACC_S = [Acc_s_X;Acc_s_Y;Acc_s_Z]';   %7.8870e-8 

Acc_n_X = sum(Acc_n.N(:,1));
Acc_n_Y = sum(Acc_n.N(:,2));
Acc_n_Z = sum(Acc_n.N(:,3));

ACC_N = [Acc_n_X;Acc_n_Y;Acc_n_Z]';   % 1.92723e-08


% 100*(0.1408- Area_tot)/0.1408;
% SRP TOT = 8.825628593774054e-08
end