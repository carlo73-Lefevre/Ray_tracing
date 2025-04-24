% In App vengono salvati gli indici riga delle facce complanari

function App = cerca_piani_app(idx_faces,center_n,App,f,Sat)

A = Sat.Vertex(Sat.Faces(:,1),:);
B = Sat.Vertex(Sat.Faces(:,2),:);
C = Sat.Vertex(Sat.Faces(:,3),:);
P = Sat.Centers_mesh;

if f < length(center_n)
    if exist('idx_faces','var')

        Fac = 1:length(center_n);

        Fi = Fac(idx_faces==0);
        Dist(1:length(center_n)) = 100;

        for  F = Fac(idx_faces==0)
            Dist(F) = Sat_piani(A(Fi(1),:),B(Fi(1),:),C(Fi(1),:),P(F,:));
        end

        idx_c  = abs(Dist) < 1e-14; % verifica quali punti sono meno distanti di una soglia
        App(idx_c) = find(idx_c ==1,1,'first');


        idh = (idx_c == 0);
        fn = sum(idx_faces);
        idx_faces(1,idx_c == 1) = true;     % mette 1 a quelli complanari


        App = cerca_piani_app(idx_faces,center_n,App,fn,Sat);

    end

end
end
