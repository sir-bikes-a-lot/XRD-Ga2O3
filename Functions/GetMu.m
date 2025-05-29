function [ mu ] = GetMu( siteAatoms, siteBatoms, n, V )
% GETMU
% Calculate the total linear absorption coefficient mu for the layer. Use
% Vegard's Law to interpolate the mass absorption coefficients.
%
% I = I0 * exp(- mu*h), where h is the layer thickness.
%
% INPUTS
% 
% siteAatoms    Atomic species and mole fractions for group III (site A).
% siteBatoms    Atomic species and mole fractions for group VI (site B).
% n             Number of molecules per unit cell; equal to 6?
% V             Volume of unit cell in Angstroms^3. Must be converted to
%               cm^3!
%
% OUTPUT is mu, the total linear absorption coefficient, in cm^-1.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Avogadro's number, N:
N = 6.022 * 10^23;

% Convert volume in Angstroms^3 to cm^3.
V = V/(10^24);

% First we need to fetch the mass absorption coefficients mu/rho for each
% element in the compound.
siteAnames = fieldnames(siteAatoms);
siteBnames = fieldnames(siteBatoms);
atomnames = [siteAnames; siteBnames];

%Initialize mu_roh and atomic mass A before entering the loop.
mu_rho = zeros(1, length(atomnames));
A = mu_rho;
mu_a = mu_rho;

% Fetch mu/rho and atomic mass A for each element. Calculate atomic 
% absorption coefficients mu_a.
for j=1:length(atomnames)
    element(j) = elements(cell2mat(atomnames(j)));
    mu_rho(j) = element(j).mu;
    A(j) = element(j).A;
    mu_a(j) = mu_rho(j)*A(j)/N;
end

% Compound is a quaternary of the form A(x)B(y)C(1-x-y)D.
% Binary1 is CD, binary2 is BD, binary3 is AD.
% Mole fractions x and y are given by:
x = siteAatoms.(cell2mat(siteAnames(1)));
y = siteAatoms.(cell2mat(siteAnames(2)));

% Apply Vegard's Law to the group III atomic absorption
% coefficients.
mu = (1 - x - y)*mu_a(3) + y*mu_a(2) + x*mu_a(1);

% Now add the group V atomic absorption coefficient.
mu = mu + mu_a(4);
mu = (n/V)*mu;   

end

