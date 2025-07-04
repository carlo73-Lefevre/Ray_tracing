function y = r2r (x)

% revolutions to radians function
   
% input

%  x = argument (revolutions; 0 <= x <= 1)

% output

%  y = equivalent x (radians; 0 <= y <= 2 pi)

% Celestial Computing with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = 2 * pi * (x - fix(x));
