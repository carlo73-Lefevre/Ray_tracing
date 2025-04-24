function Sat = centri_mesh(Sat)

% % determina i centri delle mesh
% for ii =1:length(Sat.Faces)  
%      C(ii,1) = (Sat.Vertex(Sat.Faces(ii,1),1) + Sat.Vertex(Sat.Faces(ii,2),1) +Sat.Vertex(Sat.Faces(ii,3),1))/3;
%      C(ii,2) = (Sat.Vertex(Sat.Faces(ii,1),2) + Sat.Vertex(Sat.Faces(ii,2),2) +Sat.Vertex(Sat.Faces(ii,3),2))/3;
%      C(ii,3) = (Sat.Vertex(Sat.Faces(ii,1),3) + Sat.Vertex(Sat.Faces(ii,2),3) +Sat.Vertex(Sat.Faces(ii,3),3))/3;
% end
% 
% Sat.Centri_mesh = C;

    % Calcola i centri delle mesh senza utilizzare un ciclo for
    x = Sat.Vertex(Sat.Faces(:, 1), 1) + Sat.Vertex(Sat.Faces(:, 2), 1) + Sat.Vertex(Sat.Faces(:, 3), 1);
    y = Sat.Vertex(Sat.Faces(:, 1), 2) + Sat.Vertex(Sat.Faces(:, 2), 2) + Sat.Vertex(Sat.Faces(:, 3), 2);
    z = Sat.Vertex(Sat.Faces(:, 1), 3) + Sat.Vertex(Sat.Faces(:, 2), 3) + Sat.Vertex(Sat.Faces(:, 3), 3);

    % Calcola le coordinate medie dei vertici per ottenere i centri delle facce
    C = [x, y, z] / 3;

    % Assegna i centri delle mesh all'oggetto Sat
    Sat.Centers_mesh = C;

end