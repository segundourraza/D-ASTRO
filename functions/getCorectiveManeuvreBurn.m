function [DeltaV] = getCorectiveManeuvreBurn(a, e, i, a_target, i_target, mu, r)

theta = abs(i_target - i); % Orbital inclination change
ra1 = a.*(1+e); % Apogee radius of post aerocapture orbit
at = 0.5.*(ra1+a_target); % Semi-major axis of transfer orbit
v1a = sqrt(mu.*(2./ra1-1./a)); % Velocity at apogee of post aerocapture orbit

DV_hoffman = abs(sqrt(mu) .*(sqrt(1./a_target)+sqrt(2./ra1-1./at)-sqrt(2./a_target-1./at)-sqrt(2./ra1-1./a))); % Delta V for hoffman transfer
DV_plane = 2*v1a.*sin(theta/2); % Plane change Delta V
DeltaV = (DV_plane + DV_hoffman)'; % Total Delta V

if a < 0
    vhp = sqrt(mu.*(2./r - 1./a));
    v2p = sqrt(mu.*(2./(r) - 1./a_target));
    DeltaV = abs(vhp - v2p);
end