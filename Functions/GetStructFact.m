function [ F ] = GetStructFact( siteAatoms, siteBatoms, thetaB, lambda, plane )
% GETSTRUCTFACT
% Calculate x-ray structure factor for compound. Ternary and quaternary
% compounds are handled by first calculating the atomic scattering factors
% for the constituent binaries, then using Vegard's law to interpolate
% between them.
%
% INPUTS
% siteAatoms    Structure containing the group III atoms & mole fractions.
% siteBatoms    Structure containing the group V atoms & mole fractions.
% thetaB        Bragg angle, in radians
% lambda        Wavelength of x-ray radiation, in Angstroms.
% plane         Diffraction plane, a 1-D vector of the form (h,k,l).
%
% OUTPUT is the x-ray structure factor F, a 1-D vector dependent on theta.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First we need to calculate the atomic scattering factors for each atom in
% the compound. Begin by finding f0 from the Cromer-Mann coefficients.
siteAnames = fieldnames(siteAatoms);
siteBnames = fieldnames(siteBatoms);
atomnames = [siteAnames; siteBnames];

%Initialize some variables before entering the loop.
f0 = zeros(1,length(atomnames));
f = f0;

% Fetch the Cromer-Mann coefficients from the database of elements. Then,
% use Cromer-Mann formula to calculate f0.
for j=1:length(atomnames)
    element(j) = elements(cell2mat(atomnames(j)));
    coeffs(j) = element(j).CromerMann;
    f0(j) = CromerMann(coeffs(j), thetaB, lambda);
    
    % Now add the dispersion correction factors to get atomic scattering
    % factors, f.
    f(j) = f0(j) + element(j).fprime + 1i*element(j).f2prime;
    
    % DEBUG: Finally, apply the Debye-Waller factor to the atomic 
    % scattering factors to account for mean atomic displacement due to 
    % temperature & zero-point energy.
    % f(:,j) = f(:,j).*transpose(exp(-(sin(theta)/lambda).^2 .* element(j).B));
end

% Compound is a quaternary of the form A(x)B(y)C(1-x-y)D.
% Binary1 is CD, binary2 is BD, binary3 is AD.
% Mole fractions x and y are given by:
x = siteAatoms.(cell2mat(siteAnames(1)));
y = siteAatoms.(cell2mat(siteAnames(2)));

% Apply Vegard's Law to get the quaternary atomic scattering factors.
fa = x*f(1) + y*f(2) + (1-x-y)*f(3);

% Calculate the x-ray structure factor.
F = CalcStructFact(fa, f(4), plane);

end

