function [I, I0, lamella] = Xray(sample, omega, orient, plane, offcut, delpsi, ScanType, Twotheta, lambda, thetaM, bounces, N_sigma, N_alpha, N_beta, DistAlpha, DistBeta)
% Xray
% Simulate the diffracted intensity vs incident X-ray angle (omega).
% Inputs:
% sample        Array of structures containing information about each layer
%               (lamella)
% omega         Incident x-ray angle to substrate, a 1-D vector (radians)
% orient        Miller indices (u,v,w) for substrate normal
% plane         Miller indices (h,k,l) for diffraction plane
% offcut        Substrate offcut angle w.r.t. (u,v,w), in radians
% delpsi        Azimuthal rotation w.r.t. offcut direction, in radians
% ScanType      Type of scan, either 'Omega-2Theta' or 'Rocking'
% Twotheta      If rocking curve scan, use 2theta as detector angle. (rad)
% lambda    	X-ray wavelength, for Cu_kalpha: 1.540594 Angstroms
% thetaM        Monochromator angle, in radians
% bounces       Number of monochromator bounces, usually 2 or 4
% N_sigma       Number of std. deviations for Gaussian distributions
% N_alpha       Number of crystallites in tilt space
% N_beta        Number of crystallites in azimuth space
% DistAlpha     Type of distribution in tilt space, either "Gaussian" or
%               "Uniform"
% DistBeta      Type of distribution in azimuth space, either "Gaussian" or
%               "Uniform"
%
% Outputs:
% I             Simulated X-ray intensity from sample structure
% I0            Simulated X-ray intensity from substrate only
% lamella       Array of structures for processed sample layers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pull in the binary compound and elements libraries.
names = sheetnames('binaries.xlsx');
BinariesLib = {};

for i=1:length(names)
    BinariesLib.(names(i)) = binaries(names(i));
end

names = sheetnames('elements.xlsx');
ElementsLib = {};

for i=1:length(names)
    ElementsLib.(names(i)) = elements(names(i));
end

% Loop through each layer of the sample and set up the structure for
% ProcessLamella().
N = length(sample);

a0 = Vegard(sample(1).siteAatoms, sample(1).siteBatoms, BinariesLib, 'a');          % Angstroms
b0 = Vegard(sample(1).siteAatoms, sample(1).siteBatoms, BinariesLib, 'b');          % Angstroms
c0 = Vegard(sample(1).siteAatoms, sample(1).siteBatoms, BinariesLib, 'c');          % Angstroms
beta0 = Vegard(sample(1).siteAatoms, sample(1).siteBatoms, BinariesLib, 'beta');    % degrees

% Calculate tilt between substrate normal and lattice planes. For a
% symmetric scan, tilt = 0. Otherwise scan is asymmetric.
phi = CalcTilt(orient, plane, offcut);      % in radians

% Loop through each layer of the sample and process the variables needed to
% calculate the deviation parameters.
for n=1:N
    lamella(n) = ProcessLamella(sample(n), BinariesLib, ElementsLib, omega, phi, delpsi, plane, ScanType, Twotheta, lambda, a0, b0, c0, beta0, orient, thetaM);
end

% Set up the crystallite tilt (alpha) and azimuth (beta) deviation angles
% for use in calculating deviation parameter (eta).
% Additionally, set up the weighting distributions W_alpha and W_beta.
ang_alpha = zeros(N,N_alpha);
ang_beta = zeros(N,N_beta);
W_alpha = ang_alpha;
W_beta = ang_beta;

for n=1:N
    sigma_alpha = sample(n).sigma_alpha;
    sigma_beta = sample(n).sigma_beta;
    
    switch DistAlpha
        case 'Gaussian'
            ang_alpha(n,:) = N_sigma*sigma_alpha*linspace(-1,1,N_alpha);
        case 'Uniform'
            ang_alpha(n,:) = sigma_alpha*linspace(-1,1,N_alpha);
        otherwise
            error('Invalid distribution specified for crystallite tilt, must be "Gaussian" or "Uniform"');
    end
    
    switch DistBeta
        case 'Gaussian'
            ang_beta(n,:) = N_sigma*sigma_beta*linspace(-1,1,N_beta);
        case 'Uniform'
            ang_beta(n,:) = sigma_beta*linspace(-1,1,N_beta);
        otherwise
            error('Invalid distribution specified for crystallite azimuth, must be "Gaussian" or "Uniform"');
    end

    W_alpha(n,:) = WeightDist(ang_alpha(n,:), sigma_alpha, N_alpha, DistAlpha);
    W_beta(n,:) = WeightDist(ang_beta(n,:), sigma_beta, N_beta, DistBeta);
end

% First calculate deviation parameter eta for the 1st layer (substrate).
% Calculate solution to Takagi-Taupin equation X using the Darwin-Prins
% equation (assumes infinite thickness, no dislocations).
eta = zeros(N, N_alpha, N_beta, length(omega));
X = eta;

eta0 = CalcEta(lamella(1), 0, 0);
for i=1:N_alpha
    for j=1:N_beta
        eta(1,i,j,:) = eta0;
        X(1,i,j,:) = eta0 - sign(real(eta0)).*sqrt(eta0.^2 - 1);
    end
end

% Calculate unbroadened intensity of substrate.
I0 = abs(squeeze(X(1,1,1,:))).^2;

% Now loop through and calculate eta and diffracted signal X for all
% lamellae n (except substrate) and all mosaic crystallites i,j.
% Apply mass attenuation to the diffracted signal amplitude. Assume
% mass attenuation does not influence signal phase.

for n=2:N
    for i=1:N_alpha
        for j=1:N_beta
            eta(n,i,j,:) = CalcEta(lamella(n), ang_alpha(n,i), ang_beta(n,j));
%             X(n,i,j,:) = CalcX(eta(n,i,j,:), X(n-1,i,j,:), lamella(n).T);
            X(n,i,j,:) = CalcX(squeeze(eta(n,i,j,:)), squeeze(X(n-1,i,j,:)), transpose(lamella(n).T));
            h = (10^-7)*sample(n).thickness;    % cm
            X(n,i,j,:) = squeeze(X(n,i,j,:)).*transpose(exp(-lamella(n).mu*h./(2*sin(lamella(n).Theta))));
        end
    end
end

% Calculate the total diffracted intensity from the Nth layer by summing
% over all the mosaic crystallites with the weighting functions.
I = zeros(length(omega),1);

for i=1:N_alpha
    for j=1:N_beta
        I = I + (abs(squeeze(X(N,i,j,:))).^2 * W_alpha(N,i) * W_beta(N,j));
    end
end

% Specify instrumental broadening full width half maximum (FWHM) in
% radians. Use the formula from X'Pert Epitaxy help contents, "Default 
% Convolution Conditions for Simulations". Use the angular divergence and
% wavelength spread values from the X'Pert Epitaxy sample output files.
% Assume 4-bounce Ge (220) hybrid monochromator & CuKalpha radiation.
% Monochromator angle thetaM = 22.6488 degrees for Ge (220).
thetaC = (max(omega) + min(omega))/2;           % Radians
dtheta = (pi/180) * 0.0026667 * (4/bounces);	% Radians
dlambda = 0.0004183 * (4/bounces);              % Unitless

FWHM = (2/2.64)*sqrt(dtheta^2 + dlambda^2*(tan(thetaC) - tan(thetaM)).^2);

I = InstrBroaden(I, omega, FWHM);
I0 = InstrBroaden(I0, omega, FWHM);

end

