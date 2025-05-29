function [ strained, strains ] = GetStrainedLattice(relaxed, elastic, substrate, orient, r)
% GETBIAXSTRAIN
% Calculate the out-of-plane lattice constant for the monoclinic beta-Ga2O3
% unit cell, assuming coherent strain to the substrate.
%
% INPUTS:
% relaxed       a, b, c, and beta for fully relaxed layer
% elastic       elastic constants for layer, in GPa
% substrate     a0, b0, c0, and beta0 for substrate
% orient        substrate orientation (h, k, l), a string
% r             layer relaxation, 0 <= r <= 1
%
% OUTPUTS:
% strained      ac, bc, cc, and betac for the coherently strained layer
% strains       esp1, eps2, eps3, eps5 for the coherently strianed layer
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Unpack lattice constants and elastic constants.
a = relaxed(1);
b = relaxed(2);
c = relaxed(3);
beta = relaxed(4);

a0 = substrate(1);
b0 = substrate(2);
c0 = substrate(3);
beta0 = substrate(4);

C11 = elastic(1);
C12 = elastic(2);
C13 = elastic(3);
C15 = elastic(4);
C22 = elastic(5);
C23 = elastic(6);
C25 = elastic(7);
C33 = elastic(8);
C35 = elastic(9);
C55 = elastic(10);

% Consider three substrate orientations: (010), (001), and (100).
orient_str = mat2str(orient);

switch(orient_str)
    case('[0 1 0]')
        % For coherently strained layers, the in-plane lattice parameters
        % ac, cc, and betac match the substrate.
        ac = a0;
        cc = c0;
        betac = beta0;
        
        % Calculate the four strains.
        eps1 = -(a - ac)/a;
        eps3 = cc*sin(betac*pi/180)/(c*sin(beta*pi/180)) - 1;
        eps5 = cc*cos(betac*pi/180)/(c*sin(beta*pi/180)) - ac*cos(beta*pi/180)/(a*sin(beta*pi/180));
        eps2 = -(C12*eps1 + C23*eps3 + C25*eps5)/C22;
        
        % Calculate the out-of-plane lattice parameter bc.
        bc = b*(1 + eps2);
        
    case('[0 0 1]')
        % For coherently strained layers, the in-plane lattice parameters
        % ac and bc match the substrate.
        ac = a0;
        bc = b0;
        
        % Assume that betac = beta0, justified by the small mismatch
        % between Ga2)3, Al2O3, and In2O3.
        betac = beta0;
        
        % Calculate two of the four strains.
        eps1 = -(a - ac)/a;
        eps2 = -(b - bc)/b;
        
        % Calculate the out-of-plane lattice parameter cc.
        cc = c*(C33*sin(beta*pi/180) - (C13*eps1 + C23*eps2)*sin(beta*pi/180) + (ac/a)*C35*cos(beta*pi/180));
        cc = cc/(C33*sin(betac*pi/180) + C35*cos(betac*pi/180));
        
        % Calculate the other two strains. For information only
        eps5 = cc*cos(betac*pi/180)/(c*sin(beta*pi/180)) - ac*cos(beta*pi/180)/(a*sin(beta*pi/180));
        eps3 = -(C13*eps1 + C23*eps2 + C35*eps5)/C33;
        
    case('[1 0 0]')
        % TO-DO: Need to calculate the transformed elastic constants for
        % this orientation!
        C11 = 258; C12 = 118; C13 = 139; C15 = -23;
        
        % For coherently strained layers, the in-plane lattice parameters
        % bc and cc match the substrate.
        bc = b0;
        cc = c0;
        
        % Assume that betac = beta0, justified by the small mismatch
        % between Ga2)3, Al2O3, and In2O3.
        betac = beta0;
        
        % Calculate two of the four strains.
        eps2 = -(b - bc)/b;
        eps3 = -(c - cc)/c;
        
        % Calculate the out-of-plane lattice parameter ac.
        ac = a*((C13*eps3 + C12*eps2)*sin(beta*pi/180) + (cc/c)*C15*cos(betac*pi/180) - C11*sin(beta*pi/180));
        ac = ac/(C15*cos(beta*pi/180) - C11*sin(betac*pi/180));
        
        % Calculate the other two strains. For information only
        eps5 = cc*cos(betac*pi/180)/(c*sin(beta*pi/180)) - ac*cos(beta*pi/180)/(a*sin(beta*pi/180));
        eps1 = -(C13*eps3 + C12*eps2 + C15*eps5)/C11;
        
    otherwise
%         error('ERRCODE014: Invalid substrate orientation.');
end

% Apply relaxation. TO-DO: Is this correct?
ac = (1-r)*ac + r*a;
bc = (1-r)*bc + r*b;
cc = (1-r)*cc + r*c;
betac = (1-r)*betac + r*beta;

% Pack up the coherently strained lattice parameters.
strained = [ac, bc, cc, betac];
strains = [eps1, eps2, eps3, eps5];

end

