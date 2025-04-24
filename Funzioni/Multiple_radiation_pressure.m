function F = radiation_pressure(n, Irradiance, Reflectance, Area)
% Calcola la forza di pressione di radiazione sulla n-esima superficie
%
% Inputs:
% n - numero di riflessioni
% Irradiance - Irraggiamento incidente sulla prima superficie [W/m^2]
% Reflectance - array di coefficienti di riflessione delle superfici [0-1]
% Area - area della superficie [m^2]
%
% Output:
% F - forza di pressione di radiazione [N]

% Verifica che n sia minore o uguale al numero di riflessioni precedenti
if n > length(Reflectance)
    error('Il numero di riflessioni deve essere inferiore o uguale al numero di coefficienti di riflessione');
end

% Calcola l'irraggiamento incidente sulla n-esima superficie
Irradiance_n = Irradiance * prod(Reflectance(1:n-1));

% Calcola la forza di pressione di radiazione sulla n-esima superficie
F = Irradiance_n * Area;
