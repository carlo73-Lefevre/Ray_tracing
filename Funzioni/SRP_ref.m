function  [ACC_S,ACC_N,Acc_n_tot,Acc_SRP] = SRP_ref(Sat,ray1,ray2,sun_dir,m_sat,FI_mod)

FI = 1360.8;        % ;  % W/m^2
c_const = physconst('LightSpeed');  % [m/s]

Area_tot = [];

for  ii = size(ray2,2) :-1:1  % cicla per ogni tempo

    % cicla per ogni faccia
    for ff = 1 : length(ray2(1,ii).face_source)

        face_from         = ray2(1,ii).face_from(ff);  % faccia da cui viene il raggio
        mat_from          = Sat.materials_idx(Sat.Faces(face_from,1));
        rho_from          = Sat.optical_par(mat_from).rho;

% if isempty(ray1) == 0
%     face_from         = ray1(1,ii).face_from(ff);  % faccia da cui viene il raggio
%     mat_from          = Sat.materials_idx(Sat.Faces(face_from,1));
%     rho_from = rho_from *  Sat.optical_par(mat_from).rho;
% end
        a  = Sat.Vertex(Sat.Faces(ray2(1,ii).face_source(ff),1),:);
        b  = Sat.Vertex(Sat.Faces(ray2(1,ii).face_source(ff),2),:);
        c  = Sat.Vertex(Sat.Faces(ray2(1,ii).face_source(ff),3),:);

        af  = Sat.Vertex(Sat.Faces(ray2(1,ii).face_from(ff),1),:);
        bf  = Sat.Vertex(Sat.Faces(ray2(1,ii).face_from(ff),2),:);
        cf  = Sat.Vertex(Sat.Faces(ray2(1,ii).face_from(ff),3),:);


        n       = Sat.Normals_mesh(ray2(1,ii).face_source(ff),:);
        Area_source    = calc_area(a,b,c);
        Area_from      = calc_area(af,bf,cf);
        mat_idx = Sat.materials_idx(Sat.Faces(ray2(1,ii).face_source(ff),1));
        alpha   = Sat.optical_par(mat_idx).alfa;
        rho     = Sat.optical_par(mat_idx).rho;
        delta   = Sat.optical_par(mat_idx).delta;
     % disp([string(ii),string(ff)])
        cos_b   = dot(n,-ray2(1,ii).dir(ff,:) ); % angolo tra N e il vettore satellite_sole (o raggio incidente)


        Acc_n.S(ff)       = rho_from * (-FI*FI_mod(ii)/(m_sat*c_const)) * (alpha+delta) *Area_from * cos_b;                   % SRP sun
        Acc_n.N(ff,:)     = rho_from * (-FI*FI_mod(ii)/(m_sat*c_const)) * 2*(delta/3 +rho*cos_b)*Area_from*abs(cos_b)*n;      % SRP N (delle Aree)  *FI_mod
        Acc_n.mat(ff)     = mat_idx;
        Acc_n.A(ff)       = Area_from;
        Acc_n.cos_b(ff)   = cos_b;
        Acc_n.tho_from(ff)   = rho_from;
        Area_tot             = sum([Area_tot Area_from]);

    end

Acc_s_X = sum(Acc_n.S)*-sun_dir(ii,1);
Acc_s_Y = sum(Acc_n.S)*-sun_dir(ii,2);
Acc_s_Z = sum(Acc_n.S)*-sun_dir(ii,3);

ACC_S(ii,:)  = [Acc_s_X;Acc_s_Y;Acc_s_Z]';   %7.8870e-8 

Acc_n_X = sum(Acc_n.N(:,1));
Acc_n_Y = sum(Acc_n.N(:,2));
Acc_n_Z = sum(Acc_n.N(:,3));

ACC_N(ii,:)  = [Acc_n_X;Acc_n_Y;Acc_n_Z]';   % 1.92723e-08
Acc_n_tot(ii) = Acc_n;
Acc_SRP(ii,:) = [Acc_n_X + Acc_s_X;Acc_n_Y + Acc_s_Y;Acc_n_Z + Acc_s_Z; ];

end
end