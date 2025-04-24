% trova face complanari al punto del raggio sorgente
   Vii = Sat.Vertex(Sat.Faces(ray1.face_source(1),:),:); %vertici della faccia sorgente

   plane_ii = cross(Vii(1,:)-Vii(2,:),Vii(3,:)-Vii(2,:));
   plane_ii_n = plane_ii./norm(plane_ii);

   Vii_p = Vii(1,:) - Sat.Vertex; % vettori vertice mesh al resto del solido
  

   Plane_P_Scalar_prod =  plane_ii_n*(Vii(1,:)-Vii_p)';

   idx_log = abs(Plane_P_Scalar_prod)<1e-12;    % true se appartiene al piano

   V_app = Sat.Vertex(idx_log,:); % Vertex sulla stessa faccia della sorgente