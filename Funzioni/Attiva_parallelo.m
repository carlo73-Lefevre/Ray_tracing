
function c = Attiva_parallelo(n)
% Imposta il numero massimo di worker su 24
% if exist('c','var')
%     delete(c)
% end

c = parcluster('local');
c.NumWorkers = 24;
saveProfile(c);

% Avvia il nuovo pool parallelo con 24 worker
parpool('local', n);

% Continua con il tuo codice
% ...
