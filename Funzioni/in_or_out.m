function Dentro = in_or_out(A,B,C,P)
% A B C  vertici triangolo
% P punto da valutare se interno

Dentro = true;

ab = B-A;
bc = C-B;
ca = A-C;

abap = cross(ab,P-A);
abac = cross(ab,-ca);

bcbp = cross(bc,P-B);
bcba = cross(bc,-ab);

cacp = cross(ca,P-C);
cacb = cross(ca,-bc);

val = -1e-15;

if dot(abap,abac)<= val || dot(bcbp,bcba)<= val || dot(cacp,cacb) <= val
  
     Dentro = false;

end


end