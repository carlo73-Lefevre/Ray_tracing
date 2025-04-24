function par=read_par(file_in,par)
fid=fopen(file_in,'r');

while ~feof(fid)
    linea=fgetl(fid);
    if  ~strcmp(linea(1),'%')   % controlla che non sia un commento e ci siano nella linea elementi diversi da spazio
        [attrib, remain] = strtok(linea,'=');
        remain=remain(2:end);
        [remain,wd] = strtok(remain);
        par.(attrib)=str2num(remain);
    end
end
   