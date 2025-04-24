function mater = Material_namber(filename,path)
% path = 'C:\Users\carlo\OneDrive - INAF - Istituto Nazionale di Astrofisica\G4S\Matlab_ray\Dati Caomsol\';
% path = 'D:\OneDrive - INAF - Istituto Nazionale di Astrofisica\G4S\Matlab_ray\Dati Caomsol\';
infoFile = fullfile([path,'\Dati Caomsol\'],[filename,'.xml']);
mlStruct = xml2struct(infoFile);

c = 1;
for ii = 2:1:size(mlStruct.archive.model.material,2)
    material(c).name  = mlStruct.archive.model.material{ii}.label.Attributes.label;
    material(c).tag   = mlStruct.archive.model.material{ii}.Attributes.tag;
    material(c).delta = mlStruct.archive.model.material{ii}.propertyGroup.set{1}.Attributes.value;
    material(c).rho   = mlStruct.archive.model.material{ii}.propertyGroup.set{2}.Attributes.value;
    material(c).alfa  = mlStruct.archive.model.material{ii}.propertyGroup.set{3}.Attributes.value;
    c = c+1;
end

%%
% leggi file parametri 
opt = readtable( fullfile([path,'\Dati Caomsol\'],[filename,'.dat']));

% sostituisci i valori al posto dei nomi dei parametri ottici
for op = 1:size(material,2)
    for jj =1: size(opt,1)
        if strcmp(opt.Var1{jj},material(op).delta)
            material(op).delta = opt.Var2(jj);
        end
        if strcmp(opt.Var1{jj},material(op).rho)
            material(op).rho = opt.Var2(jj);
        end
        if strcmp(opt.Var1{jj},material(op).alfa)
            material(op).alfa = opt.Var2(jj);
        end
    end
end

mater = material;
end