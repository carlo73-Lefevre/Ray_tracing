
function Area = calc_area(a,b,c)

A1 = norm(a-b);
B2 = norm(b-c);
C3 = norm(c-a);

% Dis = sort([A1,B2,C3]);
% A = Dis(1);
% B = Dis(2);
% C = Dis(3);
% 
% 
% Area = 0.25 *( (A + (B+C)) * (C-(A-B)) * (C+(A-B)) * (A+(B-C))   )^.5;
% 

    S = (A1 + B2 + C3) / 2;
    Area = sqrt(S*(S-A1)*(S-B2)*(S-C3));  % fórmula de ERONE.


end