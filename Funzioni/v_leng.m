function mag = v_leng(v)
  sv = v.* v;       %the vector with elements
  % as square of v's elements
  dp = sum(sv);     % sum of squares -- the dot product
  mag = sqrt(dp);   % magnitude
end
