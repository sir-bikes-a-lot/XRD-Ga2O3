function [ value ] = Vegard( siteAatoms, siteBatoms, param )
% VEGARD
% Applies Vegard's Law to the parameter specified by 'param'.
%
% INPUT is the structures siteAatoms and siteBatoms, which have the element
% as the field name and mole fraction as the field value.
%
% OUTPUT is the value of the parameter after Vegard's Law is applied.
%
% This function makes use of reference tables which include lattice
% parameter and elastic constants for each binary compound.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the field names (atoms) for the compound.
siteAnames = fieldnames(siteAatoms);
siteBnames = fieldnames(siteBatoms);

% Compound is a quaternary of the form A(x)B(y)C(1-x-y)D.
% Binary1 is CD, binary2 is BD, binary3 is AD.
binary1 = cell2mat(strcat(siteAnames(3), siteBnames(1)));
binary1 = binaries(binary1);
binary1 = binary1.(param);

binary2 = cell2mat(strcat(siteAnames(2), siteBnames(1)));
binary2 = binaries(binary2);
binary2 = binary2.(param);

binary3 = cell2mat(strcat(siteAnames(1), siteBnames(1)));
binary3 = binaries(binary3);
binary3 = binary3.(param);

% Mole fractions x and y are given by:
x = siteAatoms.(cell2mat(siteAnames(1)));
y = siteAatoms.(cell2mat(siteAnames(2)));

% Apply Vegard's Law.
value = (1-x-y)*binary1 + y*binary2 + x*binary3;

end

