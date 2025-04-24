
% p0 = patch('Faces',Sat0.Faces(IDX_Facce_ill(1),:),'Vertices',Sat0.Vertex);
% p0.FaceColor = [0.5 0.5 0.5]; 

ab = B-A;
bc = C-B;
ca = A-C;
%%
% il triangolo ciano (in ombra) è oscurato dal triancolo rosso che contiene
% il punto P
figure(1)
clf

pt = 55;%:length(point_indices);
idt = Tri_ok;
idp = Pt_ok;

plot3(A(idt,1),A(idt,2),A(idt,3),'o')
hold on; grid on;
plot3(B(idt,1),B(idt,2),B(idt,3),'.')
plot3(C(idt,1),C(idt,2),C(idt,3),'*')
plot3(P(point_indices(pt),1),P(point_indices(pt),2),P(point_indices(pt),3),'+')

quiver3(A(idt,1),A(idt,2),A(idt,3),ab(idt,1),ab(idt,2),ab(idt,3)) % Ab
quiver3(B(idt,1),B(idt,2),B(idt,3),bc(idt,1),bc(idt,2),bc(idt,3)) % AB
quiver3(C(idt,1),C(idt,2),C(idt,3),ca(idt,1),ca(idt,2),ca(idt,3)) % AB



% RA2 = (A(:,1)-A(1,1)).^2 + (A(:,2)-A(1,2)).^2 + (A(:,3)-A(1,3)).^2;
% RB2 = (A(:,1)-B(1,1)).^2 + (A(:,2)-B(1,2)).^2 + (A(:,3)-B(1,3)).^2;
% RC2 = (A(:,1)-C(1,1)).^2 + (A(:,2)-C(1,2)).^2 + (A(:,3)-C(1,3)).^2;

view(-90, 90);

% plot3(A(:,1),A(:,2),A(:,3),'+')
% plot3(B(:,1),B(:,2),B(:,3),'o')
% plot3(C(:,1),C(:,2),C(:,3),'*')



figure(21);
 clf;
hold on;
% plot3(A(pt,1),A(pt,2),A(pt,3),'o')
hold on; grid on;
% plot3(B(pt,1),B(pt,2),B(pt,3),'o')
% plot3(C(pt,1),C(pt,2),C(pt,3),'o')

% plot3(P(indices_of_points_inside(pt),1),P(indices_of_points_inside(pt),2),P(indices_of_points_inside(pt),3),'*')

% Mesh completa
p0 = patch('Faces', Sat0.Faces(IDX_Facce_ill,:), 'Vertices', Sat0.Vertex);
p0.FaceColor = [0.7, 0.7, 0.7];
axis equal;
p0.FaceAlpha =0.8;
% % Facce illuminate inizialmente
% p1 = patch('Faces', Sat0.Faces(IDX_Facce_ill,:), 'Vertices', Sat0.Vertex);
% p1.FaceColor = [0.2, 0.2, 1];


% triangoli in ombra 
p2 = patch('Faces', Sat0.Faces(IDX_Facce_ill(idt),:), 'Vertices', Sat0.Vertex);
p2.FaceColor = 'c';
p2.FaceAlpha =1;

% triangoli oscurante
p4 = patch('Faces', Sat0.Faces( IDX_Facce_ill(Tri_ok),:), 'Vertices', Sat0.Vertex);
p4.FaceColor = 'r';
p4.FaceAlpha =0.7;

plot3(P(point_indices(pt),1),P(point_indices(pt),2),P(point_indices(pt),3),'+')

idf = 39;pt =4;
% Faccia in ombra
p3 = patch('Faces', Sat0.Faces(IDX_Facce_ill_vect,:), 'Vertices', Sat0.Vertex);
p3.FaceColor = 'g';
p3.FaceAlpha =1;

 plot3(P(pt,1),P(pt,2),P(pt,3),'+')    % punto in ombra
 plot3(A(465 + pt,1),A(465 + pt,2),A(465 + pt,3),'*')
 plot3(B(465 + pt,1),B(465 + pt,2),B(465 + pt,3),'*')
 plot3(C(465 + pt,1),C(465 + pt,2),C(465 + pt,3),'*')

view(-90, 45);


% plot_quiv(Sat0.Centers_mesh(IDX_Facce_ill(465 + pt),:),Sat0.Normals_mesh(IDX_Facce_ill(465 + pt),:),'r')


    % id = find(idx_a>0) % indice riga della matrice dei vertici
    %
    % figure(1)
    % clf;
    % plot3(At(1),At(2),At(3),'o')
    % hold on;grid;
    % plot3(Bt(1),Bt(2),Bt(3),'o')
    % plot3(Ct(1),Ct(2),Ct(3),'o')
    % quiver3(At(1),At(2),At(3),ab(1),ab(2),ab(3)) % ab
    % quiver3(Bt(1),Bt(2),Bt(3),bc(1),bc(2),bc(3)) % bc
    % quiver3(Ct(1),Ct(2),Ct(3),ca(1),ca(2),ca(3)) % ca
    % % plot3(P1(1),P1(2),P1(3),'+')
    %
    % plot3(P1(id(1),1),P1(id(1),2),P1(id(1),3),'r+') %P1
    % quiver3(P1(id(1),1),P1(id(1),2),P1(id(1),3),vaP1(id(1),1),vaP1(id(1),2),vaP1(id(1),3),'r') % vbP1
    % quiver3(P1(1),P1(2),P1(3),vbP1(1),vbP1(2),vbP1(3),'r') % vbP1
    % quiver3(P1(1),P1(2),P1(3),vcP1(1),vcP1(2),vcP1(3),'r') % vcP1

    % id = find(idx_a>0);
    %
    %  plot3(P2(id(1),1),P2(id(1),2),P2(id(1),3),'r+') %P1
    % quiver3(P2(1),P2(2),P2(3),vaP2(1),vaP2(2),vaP2(3),'b') % vbP2
    % quiver3(P2(1),P2(2),P2(3),vbP2(1),vbP2(2),vbP2(3),'b') % vbP2
    % quiver3(P2(1),P2(2),P2(3),vcP2(1),vcP2(2),vcP2(3),'b') % vcP2
    %
    % plot3(P3(id(1),1),P3(id(1),2),P3(id(1),3),'r+') %P1
    % quiver3(P3(1),P3(2),P3(3),vaP3(1),vaP3(2),vaP3(3),'m') % vbP3
    % quiver3(P3(1),P3(2),P3(3),vbP3(1),vbP3(2),vbP3(3),'m') % vbP3
    % quiver3(P3(1),P3(2),P3(3),vcP3(1),vcP3(2),vcP3(3),'m') % vbP3
    % legend('A','B','C','ab','bc','ca',...
    %     'P1','vap1','vbP1','vcP1',...
    %     'P2','vaP2','vbP2','vcP2',...
    %     'P3','vaP3','vbP3','vcP3')
    % view(-90,22);
    % %
    %
    % figure(21);
    % clf;
    % hold on;
    %             plot3(At(1),At(2),At(3),'o')
    %             hold on;grid;
    %             plot3(Bt(1),Bt(2),Bt(3),'o')
    %             plot3(Ct(1),Ct(2),Ct(3),'o')
    %
    % A0 = Sat0.Vertex(Facce_ill_vert(c,1),:);
    % B0 = Sat0.Vertex(Facce_ill_vert(c,2),:);
    % C0 = Sat0.Vertex(Facce_ill_vert(c,3),:);
    %
    %
    %
    %             plot3(A0(1),A0(2),A0(3),'*')
    %             plot3(B0(1),B0(2),B0(3),'*')
    %             plot3(C0(1),C0(2),C0(3),'*')
    %


    % Mesh completa
    % p0 = patch('Faces', Sat0.Faces, 'Vertices', Sat0.Vertex);
    % p0.FaceColor = [0.7, 0.7, 0.7];
    % p0.FaceAlpha = 0.2
    % axis equal;
    %
    % % Facce illuminate inizialmente
    % p1 = patch('Faces', Sat0.Faces(IDX_Facce_ill,:), 'Vertices', Sat0.Vertex);
    % p1.FaceColor = 'r';
    %
    % % Facce in ombra
    % p1 = patch('Faces', Sat0.Faces(  setdiff(IDX_Facce_ill,IDX_Facce_ill_n),:), 'Vertices', Sat0.Vertex);
    % p1.FaceColor = 'g';
    % view(-90, 45);
    % camzoom(4)
    % hold off;


    %
    %             quiver3(A(c,1),A(c,2),A(c,3),ab(c,1),ab(c,2),ab(c,3)) % Ab
    %             quiver3(B(c,1),B(c,2),B(c,3),bc(c,1),bc(c,2),bc(c,3)) % AB


