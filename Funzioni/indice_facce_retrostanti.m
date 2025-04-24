function [Dot_fronte,idx_pos_fronte] = indice_facce_retrostanti(ss,sat0)

for kk = 1:size(sat0.Centri_mesh,1)  % cicla per ogni centro mesh

    V = sat0.Centri_mesh(ss,:) - sat0.Centri_mesh(kk,:);  % vettore che unisce i centri piani
    Vn = V./norm(V);
    DD(kk) = dot(-sat0.Normals_mesh(ss,:),Vn); % proiezione raggio su congiungente sorgente-mesh


    if DD(kk) > 1e-12
        clf;
        hold on;
        p3 = patch('Faces',sat0.Faces,'Vertices',sat0.Vertex); % faccia in esame per capire se sta dietro
        p3.FaceColor = [0.7,0.7,0.7]
        plot_quiv(sat0.Centri_mesh(kk,:),V,0)    % vettore congiungente
        p1 = patch('Faces',sat0.Faces(kk,:),'Vertices',sat0.Vertex); % faccia in esame per capire se sta dietro
        p1.FaceColor = 'g';
        p2 = patch('Faces',sat0.Faces(ss,:),'Vertices',sat0.Vertex); % faccia del raggio
        p2.FaceColor = 'r';
        plot_quiv(sat0.Centri_mesh(ss,:),sat0.Normals_mesh(ss,:))    % vettore Normale faccia sorgemte
        box on;
        axis equal
    end

end

idx2    = 1:size(sat0.Faces,1);

idx_pos_fronte = idx2(DD > 1e-12);    % indici di riga delle facce lato positivo
Dot_fronte     = DD(idx_pos_fronte);  % angoli tra faccia illuinata e vongiungete faccia sorgente
end