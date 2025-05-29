function [phi] = CalcTilt(orient, plane, offcut)
% CalcTilt
% Calculate tilt between substrate normal (u,v,w) and lattice plane
% (h,k,l). Tilt is zero for a symmetric scan. Otherwise, it is an
% asymmetric scan.
%
% Inputs:
% orient        Substrate orientation (u,v,w)
% plane         Lattice plane (h,k,l)
% offcut        Offcut angle from substrate (u,v,w), in radians
%
% Outputs:
% phi           Tilt angle, in radians

% Unpack orientation and plane.
u = orient(1); v = orient(2); w = orient(3);
h = plane(1); k = plane(2); l = plane(3);

% Calculate tilt angle between these two vectors using the dot product.
phi = acos((u*h + v*k +w*l)/sqrt((u^2 + v^2 +w^2)*(h^2 + k^2 + l^2)));

% Add the substrate offcut.
phi = phi + offcut;

end

