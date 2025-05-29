function [ eta ] = CalcEta( lamella, ang_alpha, ang_beta )
% CALC ETA
% Calculates the deviation parameter as a function of theta, the
% fundamental quantity for determining the X-ray diffraction pattern.
%
% INPUTS
% ang_alpha  Diffraction angle deviation for a single mosaic crystal
%            (radians).
% ang_beta   Strain deviation for a single mosaic crystal (radians).
% lamella    A structure containing the following fields:
% lamella = struct( 'Theta', theta,...              Angle between incident
%                                                   beam and (h,k,l) (rad)
%                   'delphi', delphi,...            Epilayer tilt w.r.t to
%                                                   phi (rad)
%                   'BraggAngle', thetaB,...        Bragg angle(radians)
%                   'sigma_alpha', sigma_alpha,...  Std deviation for
%                                                   diffraction angle (rad)
%                   'sigma_beta', sigma_beta,...    Std deviation for
%                                                   strain angle (rad)
%                   'Gamma', Gamma,...              Capital gamma parameter
%                   'F0', F0,...                    Structure factor for
%                                                   (000) plane
%                   'FH', FH,...                    Structure factor for
%                                                   (hkl) plane
%                   'FHbar', FHbar,...              Structure factor for
%                                                   (-h-k-l) plane
%                   'mu', mu,...                    Total linear absorption
%                                                   coefficient (cm^-1)
%                   'T', T,...                      Thickness parameter
%                                                   (unitless)
%                   'gamma0', gamma0,...            Direction cosine of
%                                                   incident beam
%                   'gammaH', gammaH,...            Direction cosine of
%                                                   diffracted beam
%                   'C', C);                        X-ray polarization
%                                                   factor
%
% OUTPUT is the deviation parameter eta for a single mosaic crystal. This
% is a 1-D vector dependent on scan angle theta only. Unitless (?).
%
% TO-DO: How does mass absorption affect this calculation? Is it only
% relevant for the calculation of X?
% TO-DO: How does temperature affect this calculation (Debye-Waller
% coefficients)? Perhaps it can be modeled in a similar manner as the
% diffraction angle & strain spread due to dislocations.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Unpack the lamella structure.
theta = lamella.Theta;
delphi = lamella.delphi;
delpsi = lamella.delpsi;
thetaB = lamella.BraggAngle;
Gamma = lamella.Gamma;
F0 = lamella.F0;
FH = lamella.FH;
FHbar = lamella.FHbar;
gamma0 = lamella.gamma0;
gammaH = lamella.gammaH;
C = lamella.C;

% Calculate eta in two parts, numerator & denominator.
eta_num = - (gamma0/gammaH)*(transpose(theta) - thetaB - (ang_alpha + delphi)*cos(ang_beta + delpsi))*sin(2*thetaB)...
          - (1/2)*(1 - gamma0/gammaH)*Gamma*F0;
eta_den = transpose(C).*Gamma.*sqrt(transpose(abs(gamma0./gammaH)).*FH.*FHbar);
eta = eta_num./eta_den;


end

