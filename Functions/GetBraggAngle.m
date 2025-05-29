function [ thetaB ] = GetBraggAngle( a, b, c, beta, lambda, plane )
% GETBRAGGANGLE
% Calculate the Bragg angle for the specified diffraction plane.
%
% INPUT
% a         a lattice constant (Angstroms)
% a         b lattice constant (Angstroms)
% c         c lattice constant (Angstroms)
% beta      beta lattice constant (degrees)
% lambda    x-ray radiation wavelength (Angstroms)
% plane     diffraction plane, specified as a vector (h,k,l)
%
% OUTPUT is the first-order Bragg angle in radians.

h = plane(1);
k = plane(2);
l = plane(3);

% Calculate interplanar spacing using auxiliary variables.
term1 = (h/a)^2;
term2 = (k/b)^2*(sin(beta*pi/180))^2;
term3 = (l/c)^2;
term4 = -2*h*l*cos(beta*pi/180)/(a*c);
invdsq = (term1 + term2 + term3 + term4)/((sin(beta*pi/180))^2);
invd = sqrt(invdsq);

% Calculate the m-th Bragg angle.
thetaB = asin(lambda*invd/2);             % radians

end

