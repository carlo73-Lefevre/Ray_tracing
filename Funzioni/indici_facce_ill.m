function idx_ill = indici_facce_ill(rif,sat0,P_rif,facce_source_idx)

%     % prodotto scalare raggio rif e Nomali (il raggio riflesso deve
%     % andare incontro alla facci quindi la proiezioene è negativa)
DD1 = rif*sat0.Normals_mesh';

%     % verifica che stia davanti
DD2 = rif*((sat0.Centers_mesh(facce_source_idx,:) - sat0.Centers_mesh)./norm((sat0.Centers_mesh(facce_source_idx,:) - sat0.Centers_mesh)))';


idx1       = 1:size(sat0.Faces,1);
idx_ill    = idx1(and(DD1 < 1e-12 , DD2 < 1e-12)); % indici facce potenz illuminate

% if ~isempty(idx_ill)
%     figure(2)
%     clf
%     hold on;
%     % facce che formano angolo positivo con raggio
%     p = patch('Faces',sat0.Faces ,'Vertices',sat0.Vertex);
%     p.FaceColor = [0.7,0.7,0.7];
%     % faccia sorgente
%     p2 = patch('Faces',sat0.Faces(idx_ill,:) ,'Vertices',sat0.Vertex);
%     p2.FaceColor = 'r';
%     % raggio riflesso in esame
%     plot_quiv(P_rif,rif)
%     box on;
%     axis equal
%        axis equal
% end

end