function GetProductionRatesBurial(lat, lon, elev)
% GetProductionRatesBurial
% Computes scaled spallogenic and muogenic production rates at the sink
% surface. Muon parameters are evaluated at z=0 (surface); depth-dependence
% is handled analytically inside the burial integral in BurialCalcComplex.
%
% Inputs:
%   lat   — Latitude [decimal degrees, negative = S]
%   lon   — Longitude [decimal degrees, negative = W]
%   elev  — Sink surface elevation [m]
%
% Outputs (saved to SinkProductionRates.mat):
%   Ps_Be, Ps_Al           — Scaled spallogenic production [atoms/g/yr]
%   Pm1_Be, Lm1_Be         — 10Be fast muon component: production and attenuation length [g/cm2]
%   Pm2_Be, Lm2_Be         — 10Be slow muon component
%   Pm1_Al, Lm1_Al         — 26Al fast muon component
%   Pm2_Al, Lm2_Al         — 26Al slow muon component

%% ===== SLHL reference production rates =====

Pslhl_Be = 3.92;    % 10Be [atoms/g/yr], Sa scaling, Borchers et al. 2016
Pslhl_Al = 28.54;   % 26Al [atoms/g/yr]

%% ===== Scaled spallogenic production (LSD, Lifton et al. 2014) =====

LSD(lat, lon, elev, 0, 1e4, 0, 10);
LSDBe = load('LSDout.mat');
tv    = LSDBe.LSDout.tv;
dt    = [tv(1), diff(tv)];
Ps_Be = Pslhl_Be * sum(LSDBe.LSDout.Be .* dt) / tv(end);

LSD(lat, lon, elev, 0, 1e4, 0, 26);
LSDAl = load('LSDout.mat');
Ps_Al = Pslhl_Al * sum(LSDAl.LSDout.Al .* dt) / tv(end);

%% ===== Muogenic production at surface (Balco 2017, two-exponential fit) =====

gmr  = -0.03417;
dtdz = 0.0065;
h    = 1013.25 * exp((gmr/dtdz) * (log(288.15) - log(288.15 - elev*dtdz)));

mc_Be = struct('Natoms', 2.006e22, 'k_neg', 0.00191*0.704*0.1828, 'sigma0', 0.280e-30);
Beout = fit_Pmu_with_exp_1A(h, 600, 2000, 2, mc_Be, 0);
Pm1_Be = Beout.P(1);  Lm1_Be = Beout.L(1);
Pm2_Be = Beout.P(2);  Lm2_Be = Beout.L(2);

mc_Al = struct('Natoms', 1.003e22, 'k_neg', 0.296*0.6559*0.0133, 'sigma0', 3.89e-30);
Alout = fit_Pmu_with_exp_1A(h, 600, 2000, 2, mc_Al, 0);
Pm1_Al = Alout.P(1);  Lm1_Al = Alout.L(1);
Pm2_Al = Alout.P(2);  Lm2_Al = Alout.L(2);

fprintf('Ps_Be=%.4f  Ps_Al=%.4f  ratio=%.4f\n', Ps_Be, Ps_Al, Ps_Al/Ps_Be);
fprintf('Pm1_Be=%.4e (L=%.1f)  Pm2_Be=%.4e (L=%.1f)\n', Pm1_Be, Lm1_Be, Pm2_Be, Lm2_Be);
fprintf('Pm1_Al=%.4e (L=%.1f)  Pm2_Al=%.4e (L=%.1f)\n', Pm1_Al, Lm1_Al, Pm2_Al, Lm2_Al);

save('SinkProductionRates.mat', 'Ps_Be', 'Ps_Al', ...
     'Pm1_Be', 'Lm1_Be', 'Pm2_Be', 'Lm2_Be', ...
     'Pm1_Al', 'Lm1_Al', 'Pm2_Al', 'Lm2_Al');
end