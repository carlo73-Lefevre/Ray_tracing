function Fcompl= triangoli_complanari(Sat)

Fcompl_n = zeros(size(Sat.Faces,1),1);

% A, B, C, D sono matrici nx3 contenenti le coordinate dei punti
A = Sat.Vertex(Sat.Faces(:,1),:);
B = Sat.Vertex(Sat.Faces(:,2),:);
C = Sat.Vertex(Sat.Faces(:,3),:);

id = 1;
for jj = 1 : size(Sat.Faces,1)

    %punto in esame
    D = Sat.Centers_mesh(jj,:);
    if Fcompl_n(jj) ~= 0
        continue
    end

    % Vettori dai vertici dei triangoli
    AB = B - A; % Vettore da A a B
    AC = C - A; % Vettore da A a C
    AD = D - A; % Vettore da A ai punti D


    % Usando il prodotto vettoriale e scalare
    cross_product = cross(AB, AC, 2);    % Prodotto vettoriale tra AB e AC
    dot_product = dot(cross_product, AD, 2);  % Prodotto scalare con AD

    % Verifica se il determinante è vicino a zero (complanarità)
    tolerance = 1e-10;  % Tolleranza per errori numerici
    complanare = abs(dot_product) < tolerance;  % 1 se i punti sono complanari, 0 altrimenti

    if Fcompl_n(jj) == 0
        Fcompl_n(complanare) = id;
        id = id+1;
    else
        Fcompl_n(complanare) = Fcompl_n(jj);
    end

    complanare =[];

end

Fcompl = Fcompl_n;

end
