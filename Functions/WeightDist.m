function [ W ] = WeightDist(x, sigma, N, type)
% WeightDist
% Calculate the mosaic crystal weighting distribution based on specified
% distribution type.
%
% Inputs:
% x         1-D vector for distribution
% sigma     Standard deviation of distribution
% N         Number of crystallites
% type      String containing distribution type
%
% Outputs:
% W         1-D vector of distribution values (weighting function)

switch(type)
    case 'Gaussian'
        if(sigma <= 0)
            sigma = 1e-5;
        end
        W = exp(-(x.^2)./(2*sigma.^2));
        W = W./sum(W);
    case 'Uniform'
        W = 1/N;
    otherwise
        error('Invalid weighting distribution specified.');
end
end

