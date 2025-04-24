% verifica se i punti appartengono al piano

function App = Sat_piani(A,B,C,P)

PL = [A' B' C']; 

A1 = [P' B' C'];
A2 = [A' P' C'];
A3 = [A' B' P'];

App = det(A1) + det(A2)+ det(A3) - det(PL);

end