function [ I_out ] =InstrBroaden( I, omega, FWHM )
% INTSTRBROADEN
% Convolve simulated x-ray intensity with a normal distribution to model
% the effect of instrumental broadening.
%
% INPUTS
% I         Simulated (unbroadened) x-ray diffraction intensity.
% omega     Scan angle range (incident x-ray beam).
% FWHM      Full-width half-maximum of x-ray machine in units of radians.
%           Accounts for both the incident and diffracted beam optics.
%
% OUTPUT is the convolution of the diffracted x-ray intensity and the
% instrumental broadening function.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Caclulate the standard deviation of the instrumental broadening normal
% distribution.
pts_per_rad = length(omega)/(omega(length(omega)) - omega(1));
width = FWHM * pts_per_rad;
sigma = width/(2*sqrt(2*log(2)));

% Set up vector of x-values for the normal distribution. Use +/-5 standard
% deviations.
N = 10*sigma;
x = linspace(-5*sigma, 5*sigma, N);

% Calculate instrumental broadening function.
f = (1/(sqrt(2*pi)*sigma))*exp((-x.^2)./(2*sigma^2));

% Convolve instrumental broadening function with diffracted x-ray intensity
% to obtain output. Trim the ends of the output to the proper length.
I_out = conv(I, f, 'same');

end

