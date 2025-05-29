function [ f0 ] = CromerMann( coeffs, thetaB, lambda )
% CROMERMANN
% Calculates the atomic scattering factor f0 using the Cromer-Mann
% expression.
%
% INPUTS are the Cromer-Mann coefficients, the Bragg angle
% thetaB, and the x-ray radiation wavelength.
%
% OUTPUT is the atomic scattering factor f0, which is dependent on theta.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize f0 and k.
f0 = 0;
k = sin(thetaB)/lambda;

% For j = 1 to j = 4, compute the sum aj * exp(-bj * sin(theta)/lambda).
for j=1:4
    field1 = strcat('a', num2str(j));
    field2 = strcat('b', num2str(j));
    
    f0 = f0 + coeffs.(field1)*exp(-coeffs.(field2)*k^2);
end

% Add on the c coefficient at the end.
f0 = f0 + coeffs.c;

end

