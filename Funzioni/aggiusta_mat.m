clear ;
filename = 'FOC001.mat';
DAT = matfile(filename);
 ray_b = DAT.ray_b;
 Sat_b = DAT.Sat_b;
ACC   = DAT.ACC;

for beta = 1:size(ACC.beta,1)
    for alfa = 1:size(ACC.alf_x,2)
        Nx = ACC.N(beta,alfa,1);
        Ny = ACC.N(beta,alfa,2);
        Nz = ACC.N(beta,alfa,3);
        N = [Nx Ny Nz];
        Sx = ACC.sun(beta,alfa,1);
        Sy = ACC.sun(beta,alfa,2);
        Sz = ACC.sun(beta,alfa,3);
        S = [Sx Sy Sz];

        SRP = (N+S);
        Sun_mag(beta,alfa) = norm(S);
        SRP_abs(beta,alfa) = norm(SRP); %cell2mat(arrayfun(@(x) norm(SRP(x,:)),[1:length(SRP)],'UniformOutput',false));
    end

end

ACC.SRP_mag = SRP_abs;
ACC.sun_mag = Sun_mag;

save(filename,'ACC','Sat_b','ray_b','-mat')