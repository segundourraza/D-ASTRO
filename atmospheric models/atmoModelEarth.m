function [T, p, rho] = atmoModelEarth(h, density_mode)
%% INFORMATION
% This function runs trajectory simulation with currernt aerothermodynamic
% model and plots a series of plots , in addition it outputs a strcuture
% with key information
% This function intends to model the martian atmosphere as a function of
% altitude.
% Data extracted from: MarsGRAM for use for computations within other .m
% files.
% The atmospheric data currently spans from 0-150km.
%   INPUTS:     -h: Altitude (m)
%               -density_mode: Determines whether the average density is used or
%               an upper/lower bound based from Mars GRAM deviations.
%                       [1] = mean atmopshere (default)
%                       [2] = lower limit of atmopshere
%                       [3] = higher limit of atmopshere
%               -full_analysis: Option for which analysis is to be done.
%   OUTPUTS:    -T: Temperature (K)
%               -P: pressure (Pa)
%               -rho: density (Kg/m3)
% h in m

if h < 25000
    T = -131.21 + 0.00299*h + 273.15;
    p = 2.488*(T/216.6).^(-11.388)*1e3;
elseif (11000 < h)
    T = -56.46+273.15;
    p = 22.65*exp(1.73 - 0.000157*h)*1e3;
else
    T = 15.04 - 0.00649*h + 273.15;
    p = 101.29 * (T/288.08).^(5.256)*1e3;
end
rho = p./(287*T);

val =0.3;
if density_mode == 2 % low
    rho = rho*(1 - val);
elseif density_mode == 3 % high
    rho = rho*(1 + val);
end