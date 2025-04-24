function S = Read_mesh_MAT_txt(filename,path)
%% test stl import
cd([path,'\Dati Caomsol\'])
S = SAT();
mater = Material_namber(filename,path);

filename_mesh = [filename,'.txt'];

T = importdata(filename_mesh,' ',10);
V = str2double([T.textdata{5}(20:end)]);
F = str2double([T.textdata{6}(20:end)]);

Vertex =  T.data(1:V ,1:3);
Faces  =  T.data(V +1:V+F,1:3);
Mats   =  T.data(V+F +1:end,1);

% % determina i centri delle mesh
% for ii =1:length(Faces)  
%      C(ii,1) = (Vertex(Faces(ii,1),1) + Vertex(Faces(ii,2),1) +Vertex(Faces(ii,3),1))/3;
%      C(ii,2) = (Vertex(Faces(ii,1),2) + Vertex(Faces(ii,2),2) +Vertex(Faces(ii,3),2))/3;
%      C(ii,3) = (Vertex(Faces(ii,1),3) + Vertex(Faces(ii,2),3) +Vertex(Faces(ii,3),3))/3;
% end
% 

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = [" ", "["];

% Specify column names and types
opts.VariableNames = ["alfaA", "VarName2", "absorptionCoefficient", "Var4"];
opts.SelectedVariableNames = ["alfaA", "VarName2", "absorptionCoefficient"];
opts.VariableTypes = ["string", "double", "categorical", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";

% Specify variable properties
opts = setvaropts(opts, ["alfaA", "Var4"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["alfaA", "absorptionCoefficient", "Var4"], "EmptyFieldRule", "auto");

% Import the data
Param = readtable([filename,'.dat'],opts);




    S.Vertex             = Vertex;
    S.Faces              = Faces;
%     S.Centri_mesh        = C;
    S.materials_idx      = Mats;
    S.optical_par        = mater;
    S.Param              = Param;


end