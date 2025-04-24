function mater = Material_namber(filename)
clear;clc;
infoFile = (['D:\OneDrive - INAF - Istituto Nazionale di Astrofisica\G4S\Matlab_ray\Dati Caomsol\'],filename,'.xml');
mlStruct = parseXML(infoFile);

chil = mlStruct.Children(2).Children;
c = 1;
for ii = 4:2:size(chil,2)
    material(c).name  = chil(ii).Children(2).Attributes.Value;
    material(c).tag   = chil(ii).Attributes(2).Value;
    material(c).name1 = chil(ii).Children(4).Children(2).Attributes(2).Value;
    material(c).name2 = chil(ii).Children(4).Children(4).Attributes(2).Value;
    material(c).name3 = chil(ii).Children(4).Children(6).Attributes(2).Value;
    c = c+1;
end

%%
opt = readtable('SAT_FEM_optic_Params.txt');

% sostituisci i valori si Delta
for op = 1:size(material,2)
    for jj =1: size(opt,1)
        if strcmp(opt.Var1{jj},material(op).name1)
            material(op).name1 = opt.Var2(jj);
        end
        if strcmp(opt.Var1{jj},material(op).name2)
            material(op).name2 = opt.Var2(jj);
        end
        if strcmp(opt.Var1{jj},material(op).name3)
            material(op).name3 = opt.Var2(jj);
        end
    end
end

mater = material;
end