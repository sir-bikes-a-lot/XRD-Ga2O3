function [ X ] = CalcX( eta, X_last, T )
% CALCX
% Calculate the diffracted x-ray signal amplitude for a single lamella and
% a single mosaic crystallite.
%
% INPUTS
% eta       Deviation parameter (unitless)
% X_last    Diffracted x-ray amplitude from last (n-1) lamella.
% T         Thickness parameter for the current (n) lamella.
%
% OUTPUT is the diffracted x-ray intensity of the present layer, X.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First, calculate the auxiliary variables S1 and S2.
S1 = (X_last - eta + sqrt(eta.^2 - 1)).*exp(-1i*T.*sqrt(eta.^2 - 1));
S2 = (X_last - eta - sqrt(eta.^2 - 1)).*exp(1i*T.*sqrt(eta.^2 - 1));

% Now calculate diffracted signal amplitude X:
X = eta + sqrt(eta.^2 - 1).*(S1 + S2)./(S1 - S2);

end

