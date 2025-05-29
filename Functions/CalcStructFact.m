function [ F ] = CalcStructFact( fA, fB, plane )
% CALCSTRUCTFACT
% Not to be confused with GetStructFact, CalcStructFact calculates the
% x-ray structure factor for a monoclinic beta-Ga2O3 crystal with group III
% atoms at site A (0.16,0.50,0.31) and (0.09,0,0.79) and group VI atoms at 
% site B (0.50,0,0.26), (0.17,0,0.56), and (0.34,0.50,0.89).
%
% INPUTS
% fA       Atomic scattering factor for group III (site A) atom.
% fB       Atomic scattering factor for group VI (site B) atom. 
% plane    DIffraction plane of interest, vector of the form (h,k,l).
%
% OUTPUT is the binary structure factor, F.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = plane(1);
k = plane(2);
l = plane(3);

% Calculate the x-ray structure factor.
F1 = fA.*exp(-2*pi*1i*(h*0.16 + k*0.50 + l*0.31));
F2 = fA.*exp(-2*pi*1i*(h*0.09 + k*0.00 + l*0.79));
F3 = fB.*exp(-2*pi*1i*(h*0.50 + k*0.00 + l*0.26));
F4 = fB.*exp(-2*pi*1i*(h*0.17 + k*0.00 + l*0.56));
F5 = fB.*exp(-2*pi*1i*(h*0.34 + k*0.50 + l*0.89));

F = F1 + F2 + F3 + F4 + F5;

end

