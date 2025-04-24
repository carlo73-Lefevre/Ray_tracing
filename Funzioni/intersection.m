function V = intersection(A,B,C,D)

syms  s t 

eqns = [ A + (B-A)*t,C + (D-C)*s];
[s,t] = solve(eqns)



end

