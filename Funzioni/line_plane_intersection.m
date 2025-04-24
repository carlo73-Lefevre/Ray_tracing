function [I,rc] = line_plane_intersection(sun_dir, Cent_mesh, sun_dir_bis, m_sun_dir5, verbose)
%% line_plane_intersection : function to compute the intersection point
% sun_dir = direzione della linea
% Cent_mesh : punto appartenente alla linea
% sun_dir_bis : normale al piano
% sun_dir5:  punto appartenete al piano

% between the (N,u) line and the (M,n) plane of the 3D space.
%
% Author & support : nicolas.douillet (at) free.fr, 2019-2020.
%
%
% Syntax
%
% [I,rc] = line_plane_intersection(sun_dir, Cent_mesh, sun_dir_bis, sun_dir5);
% [I,rc] = line_plane_intersection(sun_dir, Cent_mesh, sun_dir_bis, sun_dir5, verbose);
%
%
% Description
%
% [I,rc] = line_plane_intersection(sun_dir, N, n, M) computes the coordinates of I,
% the intersection point between the line (sun_dir,N) and the plane (n,M).
% In the most generic case, I is a point in the 3D space, but when
% the line is stricly parallel to the plane, I is the empty set, and when
% the line is included in the plane, I is a function handle corresponding
% to the system of parametric equations of the line.
%
% [I,rc] = line_plane_intersection(sun_dir, N, n, M, verbose) displays a message in
% console when verbose is set either to logical true or real numeric 1, and
% doesn't when it is set to logical false or real numeric 0.
%
%
% Principle
%
% Based on solving Descartes plane equation :
%
% ax + by + cz + d = 0, where n = [a, b, c] is a vector normal to the plane,
%
% combined with the parametric equations system of a 3D line :
%
% x(t) = x0 + at 
% y(t) = y0 + bt
% z(t) = z0 + ct
%
% where N0 = [x0, y0, z0] is a point belonging to the line, and sun_dir = [a, b, c], a vector directing this line.
%
%
% Input arguments
%
% - u : real row or column vector double. numel(u) = 3. One director vector of the parametric line.
%
% - N : real row or column vector double. numel(N) = 3. One point belonging to the line.
%
% - n : real row or column vector double. numel(n) = 3. One normal vector to the plane.
%
% - M : real row or column vector double. numel(M) = 3. One point belonging to the plane.
%
% - verbose : logical *true (1)/false(0), to enable/disable the verbose mode.
%
%
% Output arguments
%
% - I = [xI yI zI], real row or column vector double, the intersection point.
%
% - rc : return code, integer scalar doublein the set {1,2,3}.
%        0 : void / [] intersection
%        1 : point intersection (unique).
%        2 : line intersection
%
%        rc return code is necessary to distinguish between cases where
%        (N,u) line and the (M,n) plane intersection is a single point
%        and where it is the line itself.
%
%
% Example #1 : one unique intersection point
%
% n = [1 1 1];
% M = n;
% u = [1 0 0];
% N = u; % (N,u) = (OX) axis
% [I,rc] = line_plane_intersection(u, N, n, M) % one unique intersction point expected : I = [3 0 0], rc = 1
%
%
% Example #2 : line and plane are strictly // ; no intersection
%
% n = [0 0 1];
% M = [0 0 0]; % (M,n) = (XOY) plan
% u = [1 2 0];
% N = [0 0 6];
% [I,rc] = line_plane_intersection(u, N, n, M) % line strictly // plane =>  I = [], rc = 0 expected
%
%
% Example #3 : line is included in the plane
%
% n = [1 1 1];
% M = (1/3)*[1 1 1];
% u = [1 1 -2];
% N = [0.5 0.5 0];
% [I,rc] = line_plane_intersection(u, N, n, M) % line belongs to the plane, rc = 2 expected
%% Input parsing
assert(nargin > 3,'Not enough input arguments.');
assert(nargin < 6,'Too many input arguments.');
if nargin < 5    
    verbose = true;    
else    
    assert(islogical(verbose) || isreal(verbose),'verbose must be of type either logical or real numeric.');    
end
assert(isequal(size(sun_dir),size(Cent_mesh),size(sun_dir_bis),size(m_sun_dir5)),'Inputs u, M, n, and M must have the same size.');
assert(isequal(numel(sun_dir),numel(Cent_mesh),numel(sun_dir_bis),numel(m_sun_dir5),3),'Inputs u, M, n, and M must have the same number of elements (3).');
assert(isequal(ndims(sun_dir),ndims(Cent_mesh),ndims(sun_dir_bis),ndims(m_sun_dir5)),'Inputs u, M, n, and M must have the same number of dimensions.');
%% Body
% Plane offset parameter
d = -dot(sun_dir_bis,m_sun_dir5);
% Specific cases treatment
if ~dot(sun_dir_bis,sun_dir) % n & u perpendicular vectors
    if dot(sun_dir_bis,Cent_mesh) + d == 0 % N in P => line belongs to the plane
%         if verbose
%             disp('(N,u) line belongs to the (M,n) plane. Their intersection is the whole (N,u) line.');
%         end
        I = m_sun_dir5;
        rc = 2;
    else % line // to the plane
%         if verbose
%             disp('(N,u) line is parallel to the (M,n) plane. Their intersection is the empty set.');
%         end
        I = [];
        rc = 0;
    end
else
    
    % Parametric line parameter t
    t = - (d + dot(sun_dir_bis,Cent_mesh)) / dot(sun_dir_bis,sun_dir);
    
    % Intersection coordinates
    I = Cent_mesh + sun_dir*t;
    
    rc = 1;
    
end
end % line_plane_intersection