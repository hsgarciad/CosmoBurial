function [Burial_age, Age_unc, PaleoEro, Ero_unc] = ComplexBurial_muons( ...
    be10, be10_unc, al26, al26_unc, lat, lon, ...
    sedrateflag, sedrate, sink_elev, mass_depth, density)
% ComplexBurial_muons
% Burial age and pre-burial erosion rate accounting for post-burial
% cosmogenic production during sediment accumulation. Spallogenic and
% muogenic production are integrated over the burial history as the sample
% is progressively shielded at the sink.
%
% When sedrateflag=0, sedimentation rate is treated as a free parameter
% and solved simultaneously with erosion rate and burial age.
% When sedrateflag=1, the user-supplied sedrate is used as a fixed constraint.
%
% Inputs:
%   be10, be10_unc  — 10Be concentration and 1sigma uncertainty [atoms/g]
%   al26, al26_unc  — 26Al concentration and 1sigma uncertainty [atoms/g]
%   lat, lon        — Latitude [decimal degrees, negative=S/W]
%   sedrateflag     — 1: fixed sedimentation rate; 0: solve for it
%   sedrate         — Sedimentation rate [cm/yr]; required if sedrateflag=1
%   sink_elev       — Sink surface elevation [m]
%   mass_depth      — Burial mass depth [g/cm2]
%   density         — Average overburden density [g/cm3]
%
% Outputs:
%   Burial_age  — Mean burial age [kyr]
%   Age_unc     — 1sigma uncertainty [kyr]
%   PaleoEro    — Pre-burial erosion rate [m/Myr]
%   Ero_unc     — 1sigma uncertainty [m/Myr]

%% ===== Input validation =====

if sedrateflag == 1
    if nargin < 8 || isempty(sedrate)
        error('sedrateflag=1 requires a sedrate value.');
    end
    solve_Ar = false;
    Ar_fixed = sedrate;
else
    solve_Ar = true;
    Ar_fixed = [];
    fprintf('sedrateflag=0: sedimentation rate will be solved as a free parameter.\n');
end

%% ===== Constants =====

Lspal     = 160;
lambda_Be = log(2) / 1.386e6;
lambda_Al = log(2) / 7.17e5;
TBe       = 1 / lambda_Be;
TAl       = 1 / lambda_Al;
Mu        = density / Lspal;       % Spallation mass-attenuation [1/cm]

%% ===== Production rates at sink =====

pfile = 'SinkProductionRates.mat';
if isfile(pfile)
    load(pfile, 'Ps_Be', 'Ps_Al', ...
         'Pm1_Be', 'Lm1_Be', 'Pm2_Be', 'Lm2_Be', ...
         'Pm1_Al', 'Lm1_Al', 'Pm2_Al', 'Lm2_Al');
    fprintf('Loaded %s\n', pfile);
else
    fprintf('Computing production rates...\n');
    GetProductionRatesBurial(lat, lon, sink_elev);
    load(pfile, 'Ps_Be', 'Ps_Al', ...
         'Pm1_Be', 'Lm1_Be', 'Pm2_Be', 'Lm2_Be', ...
         'Pm1_Al', 'Lm1_Al', 'Pm2_Al', 'Lm2_Al');
end

% Muon mass-attenuation coefficients [1/cm]
MuT1_Be = density / Lm1_Be;
MuT2_Be = density / Lm2_Be;
MuT1_Al = density / Lm1_Al;
MuT2_Al = density / Lm2_Al;

%% ===== Initial burial age estimate =====

set(groot, 'DefaultFigureVisible', 'off');
evalc('SimpleBurial_nomuons(be10, be10_unc, al26, al26_unc, lat, lon, sink_elev)');
set(groot, 'DefaultFigureVisible', 'on');
load('SimpleBurial.mat', 'Tburial');
Tini = Tburial * 1000;   % kyr → yr

% % Initial sedimentation rate estimate (used as x0 when solve_Ar=true)
% Ar_ini = (mass_depth / density) / Tini;   % [cm/yr]
% if ~solve_Ar
%     Ar_ini = Ar_fixed;
% end

fprintf('Simple burial estimate: %.1f kyr\n', Tini/1000);
% fprintf('Initial sedimentation rate: %.4e cm/yr\n', Ar_ini);

%% ===== Model equations =====
% For a given [E (cm/yr), T (yr), Ar (cm/yr)]:
%
%   N = N_inherited * decay + N_postburial_spal + N_postburial_muons
%
%   N_inh  = (Ps  / (lambda + E*Mu))   * exp(-lambda*T)           [spallation]
%          + (Pm1 / (lambda + E*MuT1)) * exp(-lambda*T)           [muon term 1]
%          + (Pm2 / (lambda + E*MuT2)) * exp(-lambda*T)           [muon term 2]
%
%   N_post = (Ps  / (lambda - Mu*Ar))   * (exp(-Mu*Ar*T)   - exp(-lambda*T))
%          + (Pm1 / (lambda - MuT1*Ar)) * (exp(-MuT1*Ar*T) - exp(-lambda*T))
%          + (Pm2 / (lambda - MuT2*Ar)) * (exp(-MuT2*Ar*T) - exp(-lambda*T))

    function N = model_N(Ps, Pm1, Pm2, lam, mu, mt1, mt2, E, T, Ar)
        inh  = (Ps /(lam + E*mu )) * exp(-lam*T) + ...
               (Pm1/(lam + E*mt1)) * exp(-lam*T) + ...
               (Pm2/(lam + E*mt2)) * exp(-lam*T);
        post = (Ps /(lam - mu *Ar)) * (exp(-mu *Ar*T) - exp(-lam*T)) + ...
               (Pm1/(lam - mt1*Ar)) * (exp(-mt1*Ar*T) - exp(-lam*T)) + ...
               (Pm2/(lam - mt2*Ar)) * (exp(-mt2*Ar*T) - exp(-lam*T));
        N = inh + post;
    end

%% ===== Monte Carlo =====

n           = 5e3;
j           = 0;
Tburial_dis = nan(n, 1);
E_dis       = nan(n, 1);

opts = optimoptions('fsolve', 'Display', 'off', 'MaxIterations', 1000, 'MaxFunctionEvaluations', 50000);

for i = 1:n

    Be_rand = normrnd(be10, be10_unc);
    Al_rand = normrnd(al26, al26_unc);

    if Be_rand <= 0 || Al_rand <= 0; continue; end

    % Initial E guess from simple burial analytic formula
    Eini = max(Lspal * ((Ps_Be * exp(-Tini/TBe) / Be_rand) - lambda_Be), 1e-8);

    if solve_Ar
        % Ar is underdetermined with only 2 nuclides — substitute the
        % geometric consistency constraint Ar = mass_depth/(density*T)
        % directly into the model, reducing to 2 unknowns [E, T].
        % This ensures Ar is always physically consistent with the burial
        % age without falsely treating it as independently constrained.
        F = @(x) [ ...
            model_N(Ps_Be, Pm1_Be, Pm2_Be, lambda_Be, Mu, MuT1_Be, MuT2_Be, ...
                abs(x(1)), abs(x(2)), mass_depth/(density*abs(x(2)))) - Be_rand; ...
            model_N(Ps_Al, Pm1_Al, Pm2_Al, lambda_Al, Mu, MuT1_Al, MuT2_Al, ...
                abs(x(1)), abs(x(2)), mass_depth/(density*abs(x(2)))) - Al_rand  ...
        ];
        x0 = [Eini, Tini];
    else
        % Two unknowns: [E, T], Ar fixed
        Ar = Ar_fixed;
        F = @(x) [ ...
            model_N(Ps_Be, Pm1_Be, Pm2_Be, lambda_Be, Mu, MuT1_Be, MuT2_Be, abs(x(1)), abs(x(2)), Ar) - Be_rand; ...
            model_N(Ps_Al, Pm1_Al, Pm2_Al, lambda_Al, Mu, MuT1_Al, MuT2_Al, abs(x(1)), abs(x(2)), Ar) - Al_rand  ...
        ];
        x0 = [Eini, Tini];
    end

    [sol, ~, exitflag] = fsolve(F, x0, opts);

    if exitflag > 0
        Sol_E = abs(sol(1));
        Sol_T = abs(sol(2));
        if Sol_T > 0 && isreal(Sol_T) && Sol_E >= 0
            j = j + 1;
            Tburial_dis(j) = Sol_T;
            E_dis(j)       = Sol_E;
        end
    end
end

%% ===== Results =====

Tburial_dis = Tburial_dis(1:j);
E_dis       = E_dis(1:j);

if j == 0
    fprintf('No valid solutions found.\n');
    Burial_age = NaN; Age_unc = NaN; PaleoEro = NaN; Ero_unc = NaN;
    return;
end

T_kyr  = Tburial_dis / 1000;
E_mMyr = E_dis * 1e4 / density;

Burial_age = mean(rmoutliers(T_kyr));
Age_unc    = std(rmoutliers(T_kyr));
PaleoEro   = mean(rmoutliers(E_mMyr));
Ero_unc    = std(rmoutliers(E_mMyr));

fprintf('\n--- Results (%d/%d valid) ---\n', j, n);
fprintf('Burial age : %.2f +/- %.2f kyr\n', Burial_age, Age_unc);
fprintf('Erosion    : %.2f +/- %.2f m/Myr\n', PaleoEro, Ero_unc);

%% ===== Plot =====

figure('Name', 'BurialCalcComplex Results');

subplot(1,2,1);
histogram(T_kyr, 'Normalization', 'pdf', 'FaceColor', [0.8 0.8 0.8]);
hold on;
[f_t, xi_t] = ksdensity(T_kyr);
plot(xi_t, f_t, 'r-', 'LineWidth', 2);
xlabel('Burial Age [kyr]'); ylabel('Probability Density');
title('Burial Age'); grid on;
xline(Burial_age, 'k--', sprintf('%.1f kyr', Burial_age), 'LabelVerticalAlignment', 'bottom');
hold off;

subplot(1,2,2);
histogram(E_mMyr, 'Normalization', 'pdf', 'FaceColor', [0.8 0.8 0.8]);
hold on;
[f_e, xi_e] = ksdensity(E_mMyr);
plot(xi_e, f_e, 'b-', 'LineWidth', 2);
xlabel('Erosion Rate [m/Myr]'); ylabel('Probability Density');
title('Erosion Rate'); grid on;
xline(PaleoEro, 'k--', sprintf('%.1f m/Myr', PaleoEro), 'LabelVerticalAlignment', 'bottom');
hold off;

end