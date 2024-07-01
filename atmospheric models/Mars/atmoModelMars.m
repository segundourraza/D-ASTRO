function [T, P, rho] = atmoModelMars(h, density_mode, full_analysis, models)
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

if nargin <=2
    MarsGRAMdensity = models.MarsGRAMdensity;
    rho = MarsGRAMdensity(h/1e3);
else
    if density_mode == 1
        rho= models.MarsGRAMdensity(h/1e3);
    elseif density_mode ==2
        rho= models.MarsGRAMdensity_lo(h/1e3);
    elseif density_mode == 3
        rho = models.MarsGRAMdensity_hi(h/1e3);
    end
end

if nargin == 4 && full_analysis
    T = models.MarsGRAMtemperature(h/1e3);
    P = models.MarsGRAMpressure(h/1e3);
else
    T = 0;
    P = 0;
end