%% Stl import
clear;clc;

% gm = import('Mesh_Box_wing_Comsol.txt');
filename = 'Mesh_Box_wing_Comsol.txt';

text = fileread(filename);

Index1 = strfind(text, 'Coordinates');
Index2 = strfind(text, 'Elements (triangles)');

Coord_txt = text(Index1+12:Index2-1);
delta = 77;
Coord =  [str2double(Coord_txt(1:20)) ...
          str2double(Coord_txt(21:50)) ...
          str2double(Coord_txt(51:77))];
for jj = 2:Index2/delta
    Coord_dbl = [str2double(Coord_txt(1 +delta*jj:20+delta*jj)) ...
                 str2double(Coord_txt(21+delta*jj:50+delta*jj)) ...
                 str2double(Coord_txt(51+delta*jj:77+delta*jj))];
    Coord = [Coord; Coord_dbl];
end