function Sat = unisci_bodys(Sat1,Sat2)
% prende le due mesh e le unisce
% Sat1 Box 
% Sat2 Wing

numero_vert1      = length(Sat1.Vertex);
Sat = Sat1; % Box
Sat.Vertex        = [Sat1.Vertex ;       Sat2.Vertex];
Sat.Faces         = [Sat1.Faces ;        Sat2.Faces+numero_vert1];
Sat.optical_par   = [Sat1.optical_par    Sat2.optical_par];
Sat.materials_idx = [Sat1.materials_idx+2; Sat2.materials_idx + 1 + max(Sat1.materials_idx+2)];


end