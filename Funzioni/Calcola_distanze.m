
function VV_rest_sort = Calcola_distanze(sat0,idx_face_ill,ii)
% calcola le sistanz tra tutti i centri delle mesh e la mesh in esame
V          = (sat0.Centers_mesh(ii,:) - sat0.Centers_mesh(idx_face_ill,:)./norm( sat0.Centers_mesh(ii,:) - sat0.Centers_mesh(idx_face_ill,:)));
% ordina le distanze trovate
V_dis      =  sqrt(sum(V.^2,2));
V_dis(:,2) = idx_face_ill;

% vettore distanze ordinato, prima colonna distanza seconda colonna faccia 
VV_rest_sort = sortrows(V_dis,1);
end