% Load atmopsheric models from MARSGRAM data

load('MarsGRAMdensity.mat')
load('MarsGRAMdensity_hi3sd.mat')
load('MarsGRAMdensity_lo3sd.mat')
load('MarsGRAMtemperature.mat')
load('MarsGRAMpressure.mat')
load('MarsGRAMmolarwgt.mat')

% global models

models.MarsGRAMdensity = MarsGRAMdensity;
models.MarsGRAMdensity_hi = MarsGRAMdensity_hi3sd;
models.MarsGRAMdensity_lo = MarsGRAMdensity_lo3sd;
models.MarsGRAMtemperature = MarsGRAMtemperature;
models.MarsGRAMpressure = MarsGRAMpressure;
models.MarsGRAMmolarwgt = MarsGRAMmolarwgt;

simInputs.models = models;