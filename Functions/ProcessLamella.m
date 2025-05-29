function [ lamella ] = ProcessLamella( sample, BinariesLib, ElementsLib, omega, phi, delpsi, plane, ScanType, Twotheta, lambda, a0, b0, c0, beta0, orient, thetaM )
% PROCESS LAMELLA
% Calculates the following parameters for each of the N lamellae and the
% specified XRD scan range (in radians).
% INPUT is sample structure corresponding to a single layer (lamella):
%
% sample = struct( 'siteAatoms', struct('Al', 0.1, 'Ga', 0.9,...),...
%                  'siteBatoms', struct('As', 0.8, 'Sb', 0.2,...),...
%                  'thickness', h,...
%                  'delphi', delphi,...
%                  'sigma_alpha', sigma_alpha,...
%                  'sigma_beta', sigma_beta,...
%                  'relaxation', R,...
%                  'alloytype', {'binary', 'mixedIII', 'mixedV'});
%
% Also included in the inputs are:
%   BinariesLib  Library of binary compound parameters
%   ElementsLib  Library of elemental parameters
%   omega        Incident x-ray angle to substrate, a 1-D vector (radians)
%   phi          Tilt angle between substrate (u,v,w) and lattice plane
%                (h,k,l), in radians.
%   delpsi       Azimuth angular rotation, psi - psi0, in radians.
%	plane        Diffraction plane, a 3-element vector. Example: [0, 0, 4]
%   ScanType     Type of scan, either 'Omega-2Theta' or 'Rocking'
%   Twotheta     If Rocking curve scan, use 2theta as detector angle. (rad)
%   lambda       X-ray wavelength, for Cu_kalpha: 1.540594 Angstroms
%   a0           Substrate lattice constant, in Angstroms.
%   b0           Substrate lattice constant, in Angstroms.   
%   c0           Substrate lattice constant, in Angstroms.   
%   beta0        Substrate lattice constant, in degrees.
%   orient       Substrate orientation (h, k, l), a string.
%   thetaM       Bragg angle for monochromator (if applicable), in radians.
%
% NOTE that layer thickness is specified in nm, it must be converted to
% Angstroms for calculation of the thickness parameter.
%
% OUTPUT is a lamella structure of the following form, where each field is
% either a constant, or is dependent on scan range theta only:
%
% lamella = struct( 'Theta', theta,...              Angle between incident
%                                                   beam and (h,k,l) (rad)
%                   'delphi', delphi,...            Epilayer tilt w.r.t to
%                                                   phi (rad)
%                   'delpsi', delpsi,...            Epilayer azmimuthal
%                                                   rotation (rad)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Apply Vegard's Law to find the lattice constants a, b, c, and beta, and 
% the elastic constants.
a = Vegard(sample.siteAatoms, sample.siteBatoms, BinariesLib, 'a');        % Angstroms
b = Vegard(sample.siteAatoms, sample.siteBatoms, BinariesLib, 'b');        % Angstroms
c = Vegard(sample.siteAatoms, sample.siteBatoms, BinariesLib, 'c');        % Angstroms
beta = Vegard(sample.siteAatoms, sample.siteBatoms, BinariesLib, 'beta');  % degrees
C11 = Vegard(sample.siteAatoms, sample.siteBatoms, BinariesLib, 'C11');    % GPa
C12 = Vegard(sample.siteAatoms, sample.siteBatoms, BinariesLib, 'C12');    % GPa
C13 = Vegard(sample.siteAatoms, sample.siteBatoms, BinariesLib, 'C13');    % GPa
C15 = Vegard(sample.siteAatoms, sample.siteBatoms, BinariesLib, 'C15');    % GPa
C22 = Vegard(sample.siteAatoms, sample.siteBatoms, BinariesLib, 'C22');    % GPa
C23 = Vegard(sample.siteAatoms, sample.siteBatoms, BinariesLib, 'C23');    % GPa
C25 = Vegard(sample.siteAatoms, sample.siteBatoms, BinariesLib, 'C25');    % GPa
C33 = Vegard(sample.siteAatoms, sample.siteBatoms, BinariesLib, 'C33');    % GPa
C35 = Vegard(sample.siteAatoms, sample.siteBatoms, BinariesLib, 'C35');    % GPa
C55 = Vegard(sample.siteAatoms, sample.siteBatoms, BinariesLib, 'C55');    % GPa

% Pack up the lattice constants and elastic constants.
substrate = [a0, b0, c0, beta0];
relaxed = [a, b, c, beta];
elastic = [C11, C12, C13, C15, C22, C23, C25, C33, C35, C55];

% Calculate coherently strained lattice constants.
[strained, ~] = GetStrainedLattice(relaxed, elastic, substrate, orient, sample.relaxation);

% Unpack the coherently strained lattice parameters.
ac = strained(1);
bc = strained(2);
cc = strained(3);
betac = strained(4);

% If this is an interface offset layer, manually adjust the lattice
% parameter perpendicular to the substrate by the offset.
if (strcmp(sample.layertype, 'offset'))
    if (orient(1) == 1)
        ac = ac - sample.thickness;
        sample.thickness = a0;
    elseif (orient(2) == 1)
        bc = bc - sample.thickness;
        sample.thickness = b0;
    elseif (orient(3) == 1)
        cc = cc - sample.thickness;
        sample.thickness = c0;
    end
end

% Calculate Bragg angle for the strained unit cell. Specify the m-th Bragg
% reflection.
thetaB = GetBraggAngle(ac, bc, cc, betac, lambda, plane);   % Units of radians.

% Calculate unitless Gamma parameter.
re = 2.818*10^(-5);             % Classical electron radius in Angstroms.
V = ac*bc*cc*sin(betac*pi/180);	% Unit cell volume in Angstroms.
Gamma = (re*lambda^2)/(pi*V);   % Unitless parameter.

% Calculate theta, the angle between the incident X-ray beam and the
% lattice plane (h,k,l). This is offset by the scan tilt.
delphi = sample.delphi;                         % in radians
theta = omega - phi;       % in radians

% Calculate structure factors for (000), (hkl), and (-h-k-l) planes.
F0 = GetStructFact(sample.siteAatoms, sample.siteBatoms, ElementsLib, thetaB, lambda, [0,0,0]);
FH = GetStructFact(sample.siteAatoms, sample.siteBatoms, ElementsLib, thetaB, lambda, plane);
FHbar = GetStructFact(sample.siteAatoms, sample.siteBatoms, ElementsLib, thetaB, lambda, -plane);

% Calculate total linear absorption coefficient for x-ray absorption.
mu = GetMu(sample.siteAatoms, sample.siteBatoms, ElementsLib, 6, V);   % Units of cm^-1.


% Calculate X-ray polarization factor, C. Random polarization corresponds
% to thetaM = 0.
C = (1 + (cos(2*thetaM)).^2.*(cos(2*thetaB)).^2)./(1 + (cos(2*thetaM)).^2);

% Calculate unitless thickness parameter.
h = 10 * sample.thickness;   % Units of Angstroms.
gamma0 = sin(omega);

% Direction cosine of diffracted beam depends on scan type.
if (strcmp(ScanType, 'Rocking'))
    gammaH = -sin(omega - Twotheta);
elseif (strcmp(ScanType, 'Omega-2Theta'))
    gammaH = gamma0;
else
    error('ERRCODE006: Invalid scan type! Must be either "Rocking" or "Omega-2Theta".');
end

% DEBUG: This seems to be wrong!!!
% for i=length(theta)
%     T = h*C*(pi*Gamma*sqrt(FH*FHbar))/(lambda*sqrt(abs(gamma0(i)*gammaH(i))));  % Unitless.
% end
T = 0*theta;
for i=1:length(theta)
    T(i) = h*C*(pi*Gamma*sqrt(FH*FHbar))/(lambda*sqrt(abs(gamma0(i)*gammaH(i))));  % Unitless.
end

% Setup output structure "lamella".
lamella = struct('Theta', theta,...             % Angle between incident beam and (h,k,l) (rad)
                 'delphi', delphi,...           % Epilayer tilt, relative to phi (rad)
                 'delpsi', delpsi,...           % Epilayer azimuthal rotation (rad)
                 'BraggAngle', thetaB,...       % Bragg angle(radians)
                 'sigma_alpha', sample.sigma_alpha,...      % Std deviation for diffraction angle (rad)
                 'sigma_beta', sample.sigma_beta,...        % Std deviation for strain angle (rad)
                 'Gamma', Gamma,...             % Capital gamma parameter
                 'F0', F0,...                   % Structure factor for (000) plane
                 'FH', FH,...                   % Structure factor for (hkl) plane
                 'FHbar', FHbar,...             % Structure factor for (-h-k-l) plane
                 'mu', mu,...                   % Total linear absorption coefficient (cm^-1)
                 'T', T,...                     % Thickness parameter (unitless)
                 'gamma0', gamma0,...           % Direction cosine of incident beam
                 'gammaH', gammaH,...           % Direction cosine of diffracted beam
                 'C', C);                       % X-ray polarization factor

end

