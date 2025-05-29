function [ sampleout ] = SetupSample( samplein )
% SETUP ALLOY
% Process the sample structure supplied by the SampleBuilder UI for XRD
% simulation. Expands superlattices, handles interface offset layers,
% applies compound selection rules, and specifies mole fraction(s) for input
% structure "sample". We consider only the following compound types:
%
%   AB                    Binary
%   A(x)B(1-x)C           Mixed-III ternary
%   A(x)B(1-x)C(1-y)D     Mixed-III quaternary
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize a layer counter for the sampleout array.
j = 1;

% Loop through each element of the samplein array.
for i=1:length(samplein)
    % If layer type is 'single' or 'offset', assign layer to sampleout 
    % without further action.
    if(strcmp(samplein(i).layertype, 'single') || strcmp(samplein(i).layertype, 'offset'))
        sampleout(j) = samplein(i);
        j = j + 1;
        
    % If this is the first layer of a superlattice, expand the superlattice
    % in the sampleout array.
    elseif(strcmp(samplein(i).layertype, 'SL1'))
        for n=1:samplein(i).repeats
            % Assign first SL layer.
            sampleout(j + 2*(n-1)) = samplein(i);
            
            % Assign second SL layer.
            sampleout(j + 2*(n-1) + 1) = samplein(i+1);
        end
        j = j + 2*samplein(i).repeats;
        
    % Otherwise this is the second layer of a superlattice. Doing nothing
    % since we already handled the second SL layer above.
    else
    end
end

end

