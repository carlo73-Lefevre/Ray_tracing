%% verifica tra le facce potenzialment illuminat quali sono in ombra

%% Calcola le facce che non sono in ombra.
%  La funzione aggiorna l'idice delle facce illuminati. Queste sono quelle
%  che restano dal totale.
%  la logica è: controlla che la proiezione della mesh sul piano sorgente
% contenga il punto sorgente ray. Se questo punto è il j esimo relativa
% alla j-esima faccia allora non scartare la facci (false) se invece o
% relativo a un'altra faccia allora scartalo.

function IDX_Facce_ill = escludi_facce(Sat0,MAT,nray,cont)
% cont : contatore delle facce analizzate. Se alla chiamata precedente non
% è stata trovata nessuna faccia in ombra, il cont=cont+1

IDX_Facce_ill = MAT.Dis_sort_ill;       % variabile dinamica, indici facce restanti



% proiezioni del triangolino fem sul piano sorgente
% Faccia_idx   = IDX_Facce_ill(cont);   % questo indice è lo stesso per le Facce per i ray_dir_v e ray_start_p

raggio_dir   = nray.dir(cont,:);               % direzione raggio
% raggio_start = nray.p_source(cont,:);        % punto sorgente

% posizioni vertici della mesh in esame (per verifiace se è in ombra)
idx_P    = Sat0.Faces(IDX_Facce_ill,:);

% Punti della mesh sul satellite (tra quelli illuminabili)
A        = Sat0.Vertex(idx_P(:,1),:);
B        = Sat0.Vertex(idx_P(:,2),:);
C        = Sat0.Vertex(idx_P(:,3),:);

n_dir =  repmat(raggio_dir,size(IDX_Facce_ill,1), 1);

% t = n*(M-N) / (u*n) formula che clcolala distanza lungo N-M  tra
t_F_A = sum(n_dir.*(nray.p_source-A),2) ./ sum(n_dir .* n_dir,2);
t_F_B = sum(n_dir.*(nray.p_source-B),2) ./ sum(n_dir .* n_dir,2);
t_F_C = sum(n_dir.*(nray.p_source-C),2) ./ sum(n_dir .* n_dir,2);

% intersection=N+t⋅u  formula che trova il punto di intersezione sul piano
% sorgente della proiezione dei punti di mesh
F_A = A + t_F_A.*n_dir;
F_B = B + t_F_B.*n_dir;
F_C = C + t_F_C.*n_dir;

centri_mesh = Sat0.Centers_mesh(IDX_Facce_ill,:);

Dentro_idx_A =  dentro_list_vect_002(F_A,F_B,F_C,centri_mesh);
Dentro_idx_B = [];% dentro_list_vect(Sat0,IDX_Facce_ill,F_A,F_B,F_C,F_B);
Dentro_idx_C = [];% dentro_list_vect(Sat0,IDX_Facce_ill,F_A,F_B,F_C,F_C);


Dentro_idx =[Dentro_idx_A;Dentro_idx_B;Dentro_idx_C];

% Dentro_idx =  dentro_list(Sat0_ill,n_dir,nray.p_source,IDX_Facce_ill,cont,Vertex,F_A,F_B,F_C,1e-5);
 IDX_Facce_ill( Dentro_idx) = [];

% figure(21)
% clf;
% % p0 = patch('Faces',Sat0.Faces,'Vertices',Sat0.Vertex); % facce colpite dai riflessi
% % p0.FaceColor = [0.7,0.7,0.7];
% axis equal;
% hold on
% % p1 = patch('Faces',Sat0.Faces(IDX_Facce_ill,:),'Vertices',Sat0.Vertex); % facce colpite dai riflessi
% % p1.FaceColor = [0.2,0.2,1];
% p2 = patch('Faces',Sat0.Faces( IDX_Facce_ill  ,:),'Vertices',Sat0.Vertex); % facce in ombra
% p2.FaceColor = 'r';
% view(-90, 45); 
% % escludi tutte le righe ralative alle facce in ombra



% p2 = patch('Faces',Sat0.Faces(IDX_Facce_ill,:),'Vertices',Sat0.Vertex); % facce colpite dai riflessi
% p2.FaceColor = [1,0,0];
% view(-90,0);

MAT.Dis_sort_ill = IDX_Facce_ill;

end