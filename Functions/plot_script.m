A = 2.3e7;

figure;
semilogy(omega_meas + 0.0182, I_meas, 'k-'); hold on;
semilogy(omega*180/pi, A*I, 'r-'); hold off;
axis([min(omega_meas), max(omega_meas), 5e0, 1e7]);
xlabel('Omega (deg)'); ylabel('X-ray intensity (a.u.)');
title('RA905, coupled omega/2theta scan of (002)');
legend('Measurement', 'Simulation', 'Location', 'northeast');

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
import = readmatrix('C:\Users\sschaefe\OneDrive - NREL\Characterization\XRD\MCDD model\v2 Monoclinic\Development plots\v2\AlGa2O3 on Ga2O3 100 oriented.xlsx', 'Range', 'A3:F2003');
omega = import(:,1);
I1 = import(:,2);
I2 = import(:,3);
I3 = import(:,4);
I4 = import(:,5);
I5 = import(:,6);

figure;
semilogy(omega, I1, 'k-'); hold on;
semilogy(omega, I2, 'r-');
semilogy(omega, I3, 'b-');
semilogy(omega, I4, 'g-');
semilogy(omega, I5, 'm-'); hold off;
axis([min(omega), max(omega), 1e-8, 1e0]);
xlabel('Omega (deg)'); ylabel('X-ray intensity (a.u.)');
title('100 nm (AlGa)2O3 on (100) oriented Ga2O3, (400) plane');
legend('0% Al', '5% Al', '10% Al', '15% Al', '20% Al', 'Location', 'northwest');

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
import = readmatrix('C:\Users\sschaefe\OneDrive - NREL\Characterization\XRD\MCDD model\v2 Monoclinic\Development plots\v2\AlGa2O3 on Ga2O3 010 oriented.xlsx', 'Range', 'R3:U2003');
omega = import(:,1);
I1 = import(:,2);
I2 = import(:,3);
I3 = import(:,4);

figure;
semilogy(omega, I1, 'k-'); hold on;
semilogy(omega, I2, 'r-');
semilogy(omega, I3, 'b-'); hold off;
axis([min(omega), max(omega), 1e-8, 1e0]);
xlabel('Omega (deg)'); ylabel('X-ray intensity (a.u.)');
% title('100 nm (AlGa)2O3 on (010) oriented Ga2O3, (020) plane');
title('(AlGa)2O3 / Ga2O3 superlattice on (010) oriented Ga2O3, (020) plane');
legend('Single crystal', 'Tilt std dev = 0.1 deg, azimuth = 0 deg', 'Tilt std dev = 0.1 deg, uniform azimuth', 'Location', 'northeast');