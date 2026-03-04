function [Tburial, Tburial_unc, Ero, Ero_unc] = SimpleBurial_muons(be10, be10_unc, al26, al26_unc, lat, lon, elev, mass_depth)
% SimpleBurial_muons
% Computes burial age and pre-burial erosion rate from paired 10Be/26Al
% concentrations using a Monte Carlo approach. Includes post-burial
% production by muons at the specified burial depth. Recommended for
% shallow burials (mass_depth < ~500 g/cm2). For deeper samples the muon
% correction is negligible and SimpleBurial_nomuons should be used instead.
%
% Inputs:
%   be10, be10_unc  — 10Be concentration and 1sigma uncertainty [atoms/g]
%   al26, al26_unc  — 26Al concentration and 1sigma uncertainty [atoms/g]
%   lat             — Latitude [decimal degrees, negative = S]
%   lon             — Longitude [decimal degrees, negative = W]
%   elev            — Surface elevation [m]
%   mass_depth      — Burial mass depth [g/cm2]
%
% Outputs:
%   Tburial         — Mean burial age [kyr]
%   Tburial_unc     — 1sigma uncertainty [kyr]
%   Ero             — Pre-burial erosion rate [m/Myr]
%   Ero_unc         — 1sigma uncertainty [m/Myr]

%% ===== Constants =====

Pslhl_Be  = 3.92;               % 10Be SLHL production [atoms/g/yr], Borchers et al. 2016
Pslhl_Al  = 28.54;              % 26Al SLHL production [atoms/g/yr]
z         = mass_depth;         % Burial depth [g/cm2]
lambda_Be = log(2) / 1.386e6;   % 10Be decay constant [yr^-1]
lambda_Al = log(2) / 7.17e5;    % 26Al decay constant [yr^-1]
TBe       = 1 / lambda_Be;      % 10Be mean life [yr]
TAl       = 1 / lambda_Al;      % 26Al mean life [yr]
Lspal     = 160;                % Spallation attenuation length [g/cm2]
density   = 2.65;               % Quartz density for erosion rate conversion [g/cm3]

%% ===== Scaled surface production rates (LSD, Lifton et al. 2014) =====

LSD(lat, lon, elev, 0, 1e4, 0, 10);
LSDBe = load('LSDout.mat');
tv    = LSDBe.LSDout.tv;
dt    = [tv(1), diff(tv)];
Ps_Be = Pslhl_Be * sum(LSDBe.LSDout.Be .* dt) / tv(end);

LSD(lat, lon, elev, 0, 1e4, 0, 26);
LSDAl = load('LSDout.mat');
Ps_Al = Pslhl_Al * sum(LSDAl.LSDout.Al .* dt) / tv(end);

% Muon production at burial depth (Balco 2017, two-exponential fit)
gmr  = -0.03417;
dtdz = 0.0065;
h    = 1013.25 * exp((gmr/dtdz) * (log(288.15) - log(288.15 - elev*dtdz)));

mc_Be = struct('Natoms', 2.006e22, 'k_neg', 0.00191*0.704*0.1828, 'sigma0', 0.280e-30);
Beout = fit_Pmu_with_exp_1A(h, 600, 2000, 2, mc_Be, 0);
Pz_Be = sum(Beout.P .* exp(-z ./ Beout.L));

mc_Al = struct('Natoms', 1.003e22, 'k_neg', 0.296*0.6559*0.0133, 'sigma0', 3.89e-30);
Alout = fit_Pmu_with_exp_1A(h, 600, 2000, 2, mc_Al, 0);
Pz_Al = sum(Alout.P .* exp(-z ./ Alout.L));

fprintf('Ps_Be=%.4f  Ps_Al=%.4f  ratio=%.4f\n', Ps_Be, Ps_Al, Ps_Al/Ps_Be);
fprintf('Pz_Be=%.4e  Pz_Al=%.4e\n', Pz_Be, Pz_Al);

%% ===== Monte Carlo =====
% Solves for [E (g/cm2/yr), T (yr)] in:
%   N = (Ps / (lambda + E/L)) * exp(-T/Tmean)  +  Pz * Tmean * (1 - exp(-T/Tmean))
%
% Two-step solve: first without muons to get a reliable starting point,
% then with muons. This avoids the solver finding spurious solutions.

n    = 5e3;
j    = 0;
Tburial_dis = nan(n, 1);
E_dis       = nan(n, 1);

opts = optimoptions('fsolve', 'Display', 'off', 'MaxIterations', 1000, 'MaxFunctionEvaluations', 50000);

for i = 1:n

    Be_rand = normrnd(be10, be10_unc);
    Al_rand = normrnd(al26, al26_unc);

    if Be_rand <= 0 || Al_rand <= 0; continue; end

    % Initial guesses (Granger & Muzikar 2001)
    Tini = log((Al_rand/Be_rand) / (Ps_Al/Ps_Be)) / (lambda_Be - lambda_Al);
    Eini = max(Lspal * ((Ps_Be * exp(-Tini/TBe) / Be_rand) - lambda_Be), 1e-8);

    if Tini <= 0 || ~isreal(Tini) || isnan(Tini); continue; end

    % Step 1: solve without muons for a stable starting point
    F0 = @(x) [ ...
        (Ps_Al / (lambda_Al + x(1)/Lspal)) * exp(-x(2)/TAl) - Al_rand; ...
        (Ps_Be / (lambda_Be + x(1)/Lspal)) * exp(-x(2)/TBe) - Be_rand  ...
    ];
    [sol0, ~, flag0] = fsolve(F0, [Eini, Tini], opts);
    if flag0 <= 0 || sol0(2) <= 0 || sol0(1) < 0; continue; end

    % Step 2: full solve with muons, starting from the no-muon solution
    F = @(x) [ ...
        (Ps_Al / (lambda_Al + x(1)/Lspal)) * exp(-x(2)/TAl) + Pz_Al*TAl*(1 - exp(-x(2)/TAl)) - Al_rand; ...
        (Ps_Be / (lambda_Be + x(1)/Lspal)) * exp(-x(2)/TBe) + Pz_Be*TBe*(1 - exp(-x(2)/TBe)) - Be_rand  ...
    ];
    [sol, ~, exitflag] = fsolve(F, sol0, opts);

    % Enforce physical constraint: muon-corrected age must be >= no-muon age
    if exitflag > 0 && sol(2) >= sol0(2) && isreal(sol(2)) && sol(1) >= 0
        j = j + 1;
        Tburial_dis(j) = sol(2);
        E_dis(j)       = sol(1);
    end
end

%% ===== Results =====

Tburial_dis = Tburial_dis(1:j);
E_dis       = E_dis(1:j);

if j == 0
    fprintf('No valid solutions found.\n');
    Tburial = NaN; Tburial_unc = NaN; Ero = NaN; Ero_unc = NaN;
    return;
end

T_kyr  = Tburial_dis / 1000;
E_mMyr = E_dis * (1e4 / density);

Tburial     = mean(rmoutliers(T_kyr));
Tburial_unc = std(rmoutliers(T_kyr));
Ero         = mean(rmoutliers(E_mMyr));
Ero_unc     = std(rmoutliers(E_mMyr));

fprintf('\n--- Results (%d/%d valid) ---\n', j, n);
fprintf('Burial age : %.2f +/- %.2f kyr\n', Tburial, Tburial_unc);
fprintf('Erosion    : %.2f +/- %.2f m/Myr\n', Ero, Ero_unc);

%% ===== Plot =====

figure('Name', 'SimpleBurial_muons');

subplot(1,2,1);
histogram(T_kyr, 'Normalization', 'pdf', 'FaceColor', [0.8 0.8 0.8]);
hold on;
[f_t, xi_t] = ksdensity(T_kyr);
plot(xi_t, f_t, 'r-', 'LineWidth', 2);
xlabel('Burial Age [kyr]'); ylabel('Probability Density');
title('Burial Age'); grid on;
xline(Tburial, 'k--', sprintf('%.1f kyr', Tburial), 'LabelVerticalAlignment', 'bottom');
hold off;

subplot(1,2,2);
histogram(E_mMyr, 'Normalization', 'pdf', 'FaceColor', [0.8 0.8 0.8]);
hold on;
[f_e, xi_e] = ksdensity(E_mMyr);
plot(xi_e, f_e, 'b-', 'LineWidth', 2);
xlabel('Erosion Rate [m/Myr]'); ylabel('Probability Density');
title('Erosion Rate'); grid on;
xline(Ero, 'k--', sprintf('%.1f m/Myr', Ero), 'LabelVerticalAlignment', 'bottom');
hold off;

save('SimpleBurial.mat', 'Tburial', 'Tburial_unc', 'Ero', 'Ero_unc');

end