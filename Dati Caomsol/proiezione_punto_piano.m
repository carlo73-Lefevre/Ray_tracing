function P1 = proiezione_punto_piano (n,p0,l,l0)
% n : normale piano
% p0: punto piano
% l :  dir vettore
% l0: punto vettore

t = ((p0-l0).*n) /(l.*n);

 P1 = l0+l*t;

% idx_P_imp_out = (sqrt(sum(P_imp'.^2)) > 50);
% 
% 
% P_imp(idx_P_imp_out,:) = [];

end