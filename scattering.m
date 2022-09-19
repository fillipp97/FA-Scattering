%% Parameters
A = 1;
f = 2000;
w = 2*pi*f;
c = 1500; %343;
medium = "water";
k = w/c;
d = .02;
a = .002;
rho = 1000; %420; % Kg/m^-1
Z0 = rho*c;
t=5.8;

n_phi = 1000;
phi = linspace(0,2*pi,n_phi+1); 
n_r = 1000; 
r_min = -5;
r_max = 5;
r = linspace(r_min,r_max,n_r+1);
[Phi,R] = meshgrid(phi,r);
X = R.*cos(Phi);
Y = R.*sin(Phi);


Q = (4*pi.*a*A/rho*c)*(1+1/1j*k.*a);
Q = 4*pi*A/(1j*rho*w);
p_monopole = real(1j*rho*w*Q*exp(1j*(w*t-k.*R)));
p_monopole = p_monopole / max(p_monopole, [], 'all');

new_X1 = d/2 + R.*cos(Phi);
new_X2 = -d/2 + R.*cos(Phi);
newR1 = sqrt(new_X1.^2+Y.^2);
newR2 = sqrt(new_X2.^2+Y.^2);
pa_dipole = real(1j*rho*w*Q*exp(1j*(w*t-k.*newR1)));
pb_dipole = real(1j*rho*w*Q*exp(1j*(w*t-k.*newR2)));

p_dipole = pa_dipole-pb_dipole;
p_dipole = p_dipole / max(p_dipole, [], 'all');
p_plane_wave = real(A*exp(1j*(w*t-k*X)));
p_tot = real((-k^2*a^3*A*exp(1j*(w*t-k*R))*((1/3)-cos(Phi)/2)));
p_tot = p_tot / max(p_tot, [], 'all');

p_monopole_analytic = real(- k^2*a^3*A*exp(1j*(w*t-k*R))./3);
p_monopole_analytic = p_monopole_analytic / max(p_monopole_analytic, [], 'all');

p_dipole_analytic = real((k)^2*a^3*A*exp(1j*(w*t-(k)*R))*(cos(Phi)/ 2));
p_dipole_analytic = p_dipole_analytic / max(p_dipole_analytic, [], 'all');

% The analytic solution comes from the expression of the spherical wave as
% a serie as a possible subset of solutions to the differential wave
% equation under the assumption of far field

figure(1)
tcl = tiledlayout(1,3);

nexttile
h = surf(X,Y,p_monopole_analytic);
set(h,'LineStyle','none')
xlabel('X')
ylabel('Y')
zlabel('Pressure Amplitude')
title("Monopole contribution (first order)")

nexttile
h = surf(X,Y,p_dipole_analytic);
set(h,'LineStyle','none')
xlabel('X')
ylabel('Y')
zlabel('Pressure Amplitude')
title("Dipole contribution (second order)")

nexttile
h = surf(X,Y,p_tot);
set(h,'LineStyle','none')
xlabel('X')
ylabel('Y')
zlabel('Pressure Amplitude')
title("Approximation of the scattered wave")
colormap autumn
title(tcl,"Analytic solution Spherical Harmonics ka = " + k*a + " in " + medium)




% Empirical solution carried out by using the theory and summing the
% contribution of a monopole and a theoretically modelled dipole (sum of
% two opposite phased monopoles) under the assumption of far field

figure(2)
tcl2 = tiledlayout(1,3);

nexttile
h2 = surf(X,Y,-p_monopole);
set(h2,'LineStyle','none')
xlabel('X')
ylabel('Y')
zlabel('Pressure Amplitude')
title("Empirical Monopole contribution")

nexttile
h2 = surf(X,Y,p_dipole);
set(h2,'LineStyle','none')
xlabel('X')
ylabel('Y')
zlabel('Pressure Amplitude')
title("Empirical Dipole contribution ")

nexttile
h2 = surf(X,Y,p_dipole-0.3*p_monopole);
set(h2,'LineStyle','none')
xlabel('X')
ylabel('Y')
zlabel('Pressure Amplitude')
title("Empirical Approximation of the scattered wave")
colormap autumn
title(tcl2,"Empirical solution ka = " + k*a + " in " + medium)



