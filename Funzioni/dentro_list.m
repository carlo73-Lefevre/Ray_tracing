function IDX_Facce_ill = dentro_list(Sat0,n_dir,rays_source,IDX_Facce_ill,cont,Vertex,F_A,F_B,F_C,varargin)
% P0 posizione di partenza del raggio sorgente sul piano sorgente
% 
% IDX_Facce_ill righe dei Sa0.Faces illuminate 

%righe Di Sat0.Vertex che indivani i vertici delle facce illuminate
% Facce_ill_vert    = Sat0.Faces(IDX_Facce_ill,:);

% Punti della mesh sul satellite
% A        = Sat0.Vertex(Facce_ill_vert(:,1),:);
% B        = Sat0.Vertex(Facce_ill_vert(:,2),:);
% C        = Sat0.Vertex(Facce_ill_vert(:,3),:);

% t = n*(M-N) / (u*n) formula che clcolala distanza lungo N-M  tra
% t_F_A = sum(n_dir.*(rays_source-A),2) ./ sum(n_dir .* n_dir,2);
% t_F_B = sum(n_dir.*(rays_source-B),2) ./ sum(n_dir .* n_dir,2);
% t_F_C = sum(n_dir.*(rays_source-C),2) ./ sum(n_dir .* n_dir,2);

% intersection=N+t⋅u  formula che trova il punto di intersezione sul piano
% sorgente
A = F_A;% A + t_F_A.*n_dir;
B = F_B;% B + t_F_B.*n_dir;
C = F_C;% BC + t_F_C.*n_dir;


% punto sorgente del raggio su piano sorgente

    ba = B -A ;
    bc = B -C ;
    ca = C -A ;
    Area_sect = vecnorm((cross(ba,ca))/2,2,2) ;  %calc_area(ba,bc,ca);

for c = cont:size(rays_source,1)
    ray_sourc_n_mat = repmat(rays_source(c,:),size(A,1), 1);

    Area_1    = vecnorm((cross(ba,A-ray_sourc_n_mat))/2,2,2);  %calc_area(ba,A-P,B-P);
    Area_2    = vecnorm((cross(bc,B-ray_sourc_n_mat))/2,2,2);  %calc_area(bc,B-P,C-P);
    Area_3    = vecnorm((cross(ca,C-ray_sourc_n_mat))/2,2,2);  %calc_area(ca,C-P,A-P);

    % escludi tutte le mesh che contengono questo punto treanne la prima
        riga_ray_sources_in = ((Area_1+Area_2+Area_3)- Area_sect) <= 1e-3; %


end

    [l1, ~]   = find(riga_ray_sources_in == 1);

    IDX_Facce_ill(l1(2:end))  =  [];
    rays_source(l1(2:end),:)  =  [];
    n_dir(l1(2:end),:)        =  [];
    A(l1(2:end),:)                   =  [];
    B(l1(2:end),:)                   =  [];
    C(l1(2:end),:)                   =  [];

if length(l1) == 1
    cont = cont + 1;
else 
    cont = 1;
end
    size(IDX_Facce_ill)
    
    IDX_Facce_ill =  dentro_list(Sat0,n_dir,rays_source,IDX_Facce_ill,cont,Vertex,A,B,C,varargin{1,1});



end