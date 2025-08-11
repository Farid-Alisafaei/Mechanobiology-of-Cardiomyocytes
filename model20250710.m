clear all
clc

% cardiomyocyte constants
E_t0 = 0.25;    % initial stiffness of cyto in tension kPa (cardiomycocyte)
E_c = 1;        % stiffness of cyto in compression kPa (cardiomycocyte)
beta = 32;      % CM chemical stiffness parameter
alpha = 31;     % CM Feedback feedback parameter
C0_base = 0.025;    % CM base contractility
amplitude = 0.037;  % CM contractility amplitude
n = 1;      % linear (should be 1; contractility increases linearly with stress)
TGF = 1;    % controls the contractility level of cardiomyocyte (should be 1; TGF-beta only affects the fibroblast directly)
l = 45;     % stiffening parameter of cyto in tension in cardiomycocyte
m = 0.7;    % stiffening parameter of cyto in tension in cardiomycocyte

% fibroblast constants
E_t0f = 10;     % stiffness of cyto in tension kPa (fibroblast)
E_cf = 1;       % stiffness of cyto in compression kPa (fibroblast)
betaf = 32;     % fibroblast chemical stiffness parameter
alphaf = 31;    % fibroblast Feedback feedback parameter
Rho0 = 0.01;    % fibroblast base contractility
TGFF = 1;       % controls the contractility level (effect of TGF-beta on fibroblast; 1 for untreated)
lf = 0;         % stiffening parameter of cyto in tension in fibroblast
mf = 0;         % stiffening parameter of cyto in tension in fibroblast

% collagen constants
E_col0 = 10;    % stiffness of ECM kPa
lcol = 0;       % stiffening parameter of ECM
mcol = 0.5;     % stiffening parameter of ECM
yita = 50;      % ECM viscosity

% post constants
E_p = 1000000000; % rigid posts in stretching experiment


% stability criteria check
if alpha <= (1 / E_c)
    fprintf('Stability criteria violated\n');
end

if beta <= alpha
    fprintf('Stability criteria violated\n');
end

if alphaf <= (1 / E_cf)
    fprintf('Stability criteria violated\n');
end

if betaf <= alphaf
    fprintf('Stability criteria violated\n');
end

%% poking experiment (stages 1 & 2)
% time for plot
time = 0:0.01:50;
dt = time(2) - time(1);

% initial conditions
epsilon_yita_prev = 0;

solutions = [];

% guesses for fsolve
initial_guess = zeros(1, 15);
initial_guess(8) = epsilon_yita_prev;
initial_guess(10) = C0_base;

% stretch magnitudes for plot 1
epsilon_ext_1 = 0.01;
epsilon_ext_2 = 0.3;

for i = 1:length(time)
    current_time = time(i);

    % increase epsilon_ext for stretching
    if current_time < 32
        epsilon_ext_current = epsilon_ext_1;
    elseif current_time >= 32 && current_time <= 34
        epsilon_ext_current = epsilon_ext_1 + (epsilon_ext_2 - epsilon_ext_1) * (current_time - 32) / (34 - 32);
    else
        epsilon_ext_current = epsilon_ext_2;
    end

    options = optimset('Display', 'off');
    [sol, ~] = fsolve(@(variables) system_of_equations(variables, current_time, epsilon_yita_prev, C0_base, amplitude,...
        E_p, E_t0, E_c, beta, alpha, yita, n, epsilon_ext_current, dt, l, m, TGF, E_t0f, E_cf, betaf, alphaf, Rho0, TGFF, lf, mf, E_col0, lcol, mcol), initial_guess, options);

    solutions = [solutions; sol];

    % update epsilon_yita_prev for the next time step
    epsilon_yita_prev = sol(8);

    % update the initial guess for the next iteration
    initial_guess = sol;
end


%% plot Sigma and force for poking experiment (stages 1 & 2)
figure;
plot(time, solutions(:, 5));
title('Sigma');
xlabel('Time (s)');
ylabel('Sigma (kPa)');

area = 3;   % sample area in poking experiment representative (mm2)
force = solutions(:, 5) * area;

figure;
plot(time, force, 'k-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Force (mN)');
title('Force');
grid on;

%% poking experiment (stages 1 & 2) with vs. without contractility

solutions_passive = [];

% initial conditions
epsilon_yita_prev = 0;
initial_guess = zeros(1, 15);
initial_guess(8) = epsilon_yita_prev;
initial_guess(10) = 0;  % C0_base = 0 for passive case

% passive parameters (no contractility)
C0_base_passive = 0;
amplitude_passive = 0;

for i = 1:length(time)
    current_time = time(i);

    % increase epsilon_ext for stretching
    if current_time < 32
        epsilon_ext_current = epsilon_ext_1;
    elseif current_time >= 32 && current_time <= 34
        epsilon_ext_current = epsilon_ext_1 + (epsilon_ext_2 - epsilon_ext_1) * (current_time - 32) / (34 - 32);
    else
        epsilon_ext_current = epsilon_ext_2;
    end

    options = optimset('Display', 'off');
    [sol, ~] = fsolve(@(variables) system_of_equations(variables, current_time, epsilon_yita_prev, ...
        C0_base_passive, amplitude_passive, E_p, E_t0, E_c, beta, alpha, yita, n, ...
        epsilon_ext_current, dt, l, m, TGF, E_t0f, E_cf, betaf, alphaf, ...
        Rho0, TGFF, lf, mf, E_col0, lcol, mcol), initial_guess, options);

    solutions_passive = [solutions_passive; sol];

    epsilon_yita_prev = sol(8);
    initial_guess = sol;
end

% plot sigma with vs. without contractility
figure;
plot(time, solutions(:, 5), 'r-', 'LineWidth', 2); hold on;
plot(time, solutions_passive(:, 5), 'k--', 'LineWidth', 2);
xlabel('Time');
ylabel('Stress (Sigma)');
title('Stress (with and without contractility)');
legend('with contractility', 'without contractility (C₀ = 0, amplitude = 0)');
grid on;


%% Force in poking experiment (2 stages) without strain stiffening

% time for plot
time = 0:0.01:50;
dt = time(2) - time(1);

l_nostiff = 0; % without stiffening

% initial conditions
epsilon_yita_prev = 0;

solutions = [];

% guesses for fsolve
initial_guess = zeros(1, 15);
initial_guess(8) = epsilon_yita_prev;
initial_guess(10) = C0_base;

% stretch magnitudes for plot 1
epsilon_ext_1 = 0.01;
epsilon_ext_2 = 0.3;

for i = 1:length(time)
    current_time = time(i);

    % increase epsilon_ext for stretching
    if current_time < 32
        epsilon_ext_current = epsilon_ext_1;
    elseif current_time >= 32 && current_time <= 34
        epsilon_ext_current = epsilon_ext_1 + (epsilon_ext_2 - epsilon_ext_1) * (current_time - 32) / (34 - 32);
    else
        epsilon_ext_current = epsilon_ext_2;
    end

    options = optimset('Display', 'off');
    [sol, ~] = fsolve(@(variables) system_of_equations(variables, current_time, epsilon_yita_prev, C0_base, amplitude,...
        E_p, E_t0, E_c, beta, alpha, yita, n, epsilon_ext_current, dt, l_nostiff, m, TGF, E_t0f, E_cf, betaf, alphaf, Rho0, TGFF, lf, mf, E_col0, lcol, mcol), initial_guess, options);

    solutions = [solutions; sol];

    % update epsilon_yita_prev for the next time step
    epsilon_yita_prev = sol(8);

    % update the initial guess for the next iteration
    initial_guess = sol;
end

area = 3;
force = solutions(:, 5) * area;

figure;
plot(time, force, 'k-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Force (mN)');
title('Force vs. Time without Strain Stiffening');
grid on;



%% Setup for baseline and amplitude calculations
epsilon_ext_values = linspace(0.0, 0.25, 26);   % defines the stretch levels
time_range = linspace(0, 50, 200);              % simulation time: 0–20 s
dt = time_range(2) - time_range(1);             % time step
time_steady = 45;                               % time when visocsity vanishes for calculating base and amplitude


%% baseline of Sigma in stretching experiment
sigma_min_values = zeros(size(epsilon_ext_values));
sigma_max_values = zeros(size(epsilon_ext_values));

for j = 1:length(epsilon_ext_values)
    epsilon_ext_current = epsilon_ext_values(j);
    sigma_values_over_time = [];
    time_log = [];

    % initial conditions
    epsilon_yita_prev = 0;
    initial_guess = zeros(1, 15);
    initial_guess(8) = epsilon_yita_prev;
    initial_guess(10) = C0_base;

    for t = time_range
        options = optimoptions('fsolve', 'Display', 'none');
        sol = fsolve(@(variables) system_of_equations(variables, t, epsilon_yita_prev, ...
            C0_base, amplitude, E_p, E_t0, E_c, beta, alpha, yita, n, ...
            epsilon_ext_current, dt, l, m, TGF, ...
            E_t0f, E_cf, betaf, alphaf, Rho0, TGFF, lf, mf, ...
            E_col0, lcol, mcol), initial_guess, options);

        sigma_values_over_time(end + 1) = sol(5);
        time_log(end + 1) = t;
        epsilon_yita_prev = sol(8);
        initial_guess = sol;
    end

    % only keep Sigma values after reaching steady
    steady_indices = time_log >= time_steady;
    sigma_min_values(j) = min(sigma_values_over_time(steady_indices));
    sigma_max_values(j) = max(sigma_values_over_time(steady_indices));
end

figure;
plot(epsilon_ext_values, sigma_min_values, 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Epsilon_{ext}');
ylabel('Baseline stress');
title('Baseline stress');
grid on;

figure;
plot(epsilon_ext_values, sigma_max_values, 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Epsilon_{ext}');
ylabel('Peak stress');
title('Peak stress');
grid on;

%% Baseline of Sigma in Stretching with and without contractility
sigma_min_passive = zeros(size(epsilon_ext_values));

C0_base_passive = 0;
amplitude_passive = 0;

for j = 1:length(epsilon_ext_values)
    epsilon_ext_current = epsilon_ext_values(j);
    sigma_values_passive = [];
    time_log = [];

    % initial conditions
    epsilon_yita_prev = 0;
    initial_guess = zeros(1, 15);
    initial_guess(8) = epsilon_yita_prev;
    initial_guess(10) = C0_base_passive;

    for t = time_range
        options = optimoptions('fsolve', 'Display', 'none');
        sol = fsolve(@(variables) system_of_equations(variables, t, epsilon_yita_prev, ...
            C0_base_passive, amplitude_passive, E_p, E_t0, E_c, beta, alpha, yita, n, ...
            epsilon_ext_current, dt, l, m, TGF, ...
            E_t0f, E_cf, betaf, alphaf, Rho0, TGFF, lf, mf, ...
            E_col0, lcol, mcol), initial_guess, options);

        sigma_values_passive(end + 1) = sol(5);
        time_log(end + 1) = t;
        epsilon_yita_prev = sol(8);
        initial_guess = sol;
    end

    % only steady values
    steady_indices = time_log >= time_steady;
    sigma_min_passive(j) = min(sigma_values_passive(steady_indices));
end

% plot
figure;
plot(epsilon_ext_values, sigma_min_values, 'b-o', 'LineWidth', 2, 'MarkerSize', 6); hold on;
plot(epsilon_ext_values, sigma_min_passive, 'r-o', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Epsilon_{ext}');
ylabel('Baseline Stress');
title('Baseline Stress with and without contractility');
legend('with contractility', 'without contractility (C₀ = 0, amplitude = 0)', 'Location', 'best');
grid on;



%% amplitude of Sigma in sretching experiment
sigma_max_values = zeros(size(epsilon_ext_values));
sigma_diff_values = zeros(size(epsilon_ext_values));  % amplitude = max - min

for j = 1:length(epsilon_ext_values)
    epsilon_ext_current = epsilon_ext_values(j);
    sigma_values_over_time = [];
    time_log = [];

    %  initial conditions
    epsilon_yita_prev = 0;
    initial_guess = zeros(1, 15);
    initial_guess(8) = epsilon_yita_prev;
    initial_guess(10) = C0_base;

    for t = time_range
        options = optimoptions('fsolve', 'Display', 'none');
        sol = fsolve(@(variables) system_of_equations(variables, t, epsilon_yita_prev, ...
            C0_base, amplitude, E_p, E_t0, E_c, beta, alpha, yita, n, ...
            epsilon_ext_current, dt, l, m, TGF, ...
            E_t0f, E_cf, betaf, alphaf, Rho0, TGFF, lf, mf, ...
            E_col0, lcol, mcol), initial_guess, options);

        sigma_values_over_time(end + 1) = sol(5);
        time_log(end + 1) = t;
        epsilon_yita_prev = sol(8);
        initial_guess = sol;
    end

    % only keep Sigma values after steady
    steady_indices = time_log >= time_steady;
    sigmas_steady = sigma_values_over_time(steady_indices);

    sigma_max_values(j) = max(sigmas_steady);
    sigma_min_values(j) = min(sigmas_steady);
    sigma_diff_values(j) = sigma_max_values(j) - sigma_min_values(j);
end

figure;
plot(epsilon_ext_values, sigma_diff_values, 'k-o', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Epsilon_{ext}');
ylabel('stress amplitude');
title('stress amplitude');
grid on;


%% amplitude of Sigma in sretching experiment with and without contractility
sigma_max_passive = zeros(size(epsilon_ext_values));
sigma_min_passive = zeros(size(epsilon_ext_values));
sigma_diff_passive = zeros(size(epsilon_ext_values));

for j = 1:length(epsilon_ext_values)
    epsilon_ext_current = epsilon_ext_values(j);
    sigma_values_over_time = [];
    time_log = [];

    %  initial conditions
    epsilon_yita_prev = 0;
    initial_guess = zeros(1, 15);
    initial_guess(8) = epsilon_yita_prev;
    initial_guess(10) = 0;  % passive C₀ = 0

    for t = time_range
        options = optimoptions('fsolve', 'Display', 'none');
        sol = fsolve(@(variables) system_of_equations(variables, t, epsilon_yita_prev, ...
            0, 0, E_p, E_t0, E_c, beta, alpha, yita, n, ...
            epsilon_ext_current, dt, l, m, TGF, ...
            E_t0f, E_cf, betaf, alphaf, Rho0, TGFF, lf, mf, ...
            E_col0, lcol, mcol), initial_guess, options);

        sigma_values_over_time(end + 1) = sol(5);
        time_log(end + 1) = t;
        epsilon_yita_prev = sol(8);
        initial_guess = sol;
    end

    % steady-state 
    steady_indices = time_log >= time_steady;
    sigmas_steady = sigma_values_over_time(steady_indices);

    sigma_max_passive(j) = max(sigmas_steady);
    sigma_min_passive(j) = min(sigmas_steady);
    sigma_diff_passive(j) = sigma_max_passive(j) - sigma_min_passive(j);
end

% plot 
figure;
plot(epsilon_ext_values, sigma_diff_values, 'k-o', 'LineWidth', 2, 'MarkerSize', 6); hold on;
plot(epsilon_ext_values, sigma_diff_passive, 'r-o', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Epsilon_{ext}');
ylabel('stress amplitude');
title('stress amplitude with and without contractility');
legend('with contractility', 'without contractility', 'Location', 'best');
grid on;


%% baseline of Sigma in stretching experiment without stiffening
sigma_min_values = zeros(size(epsilon_ext_values));
sigma_max_values = zeros(size(epsilon_ext_values));
l_nostiff = 0;

for j = 1:length(epsilon_ext_values)
    epsilon_ext_current = epsilon_ext_values(j);
    sigma_values_over_time = [];
    time_log = [];

    %  initial conditions
    epsilon_yita_prev = 0;
    initial_guess = zeros(1, 15);
    initial_guess(8) = epsilon_yita_prev;
    initial_guess(10) = C0_base;

    for t = time_range
        options = optimoptions('fsolve', 'Display', 'none');
        sol = fsolve(@(variables) system_of_equations(variables, t, epsilon_yita_prev, ...
            C0_base, amplitude, E_p, E_t0, E_c, beta, alpha, yita, n, ...
            epsilon_ext_current, dt, l_nostiff, m, TGF, ...
            E_t0f, E_cf, betaf, alphaf, Rho0, TGFF, lf, mf, ...
            E_col0, lcol, mcol), initial_guess, options);

        sigma_values_over_time(end + 1) = sol(5);
        time_log(end + 1) = t;
        epsilon_yita_prev = sol(8);
        initial_guess = sol;
    end

    % only keep Sigma values after t = 15 (steady)
    steady_indices = time_log >= time_steady;
    sigma_min_values(j) = min(sigma_values_over_time(steady_indices));
    sigma_max_values(j) = max(sigma_values_over_time(steady_indices));
end

figure;
plot(epsilon_ext_values, sigma_min_values, 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Epsilon_{ext}');
ylabel('Baseline stress');
title('Baseline stress');
grid on;

figure;
plot(epsilon_ext_values, sigma_max_values, 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Epsilon_{ext}');
ylabel('Max stress');
title('Max stress');
grid on;

%% amplitude of Sigma in stretching experiment without stiffening
sigma_max_values = zeros(size(epsilon_ext_values));
sigma_diff_values = zeros(size(epsilon_ext_values));  % amplitude = max - min

for j = 1:length(epsilon_ext_values)
    epsilon_ext_current = epsilon_ext_values(j);
    sigma_values_over_time = [];
    time_log = [];

    %  initial conditions
    epsilon_yita_prev = 0;
    initial_guess = zeros(1, 15);
    initial_guess(8) = epsilon_yita_prev;
    initial_guess(10) = C0_base;

    for t = time_range
        options = optimoptions('fsolve', 'Display', 'none');
        sol = fsolve(@(variables) system_of_equations(variables, t, epsilon_yita_prev, ...
            C0_base, amplitude, E_p, E_t0, E_c, beta, alpha, yita, n, ...
            epsilon_ext_current, dt, l_nostiff, m, TGF, ...
            E_t0f, E_cf, betaf, alphaf, Rho0, TGFF, lf, mf, ...
            E_col0, lcol, mcol), initial_guess, options);

        sigma_values_over_time(end + 1) = sol(5);
        time_log(end + 1) = t;
        epsilon_yita_prev = sol(8);
        initial_guess = sol;
    end

    % only keep Sigma values after steady
    steady_indices = time_log >= time_steady;
    sigmas_steady = sigma_values_over_time(steady_indices);

    sigma_max_values(j) = max(sigmas_steady);
    sigma_min_values(j) = min(sigmas_steady);
    sigma_diff_values(j) = sigma_max_values(j) - sigma_min_values(j);
end

figure;
plot(epsilon_ext_values, sigma_diff_values, 'k-o', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Epsilon_{ext}');
ylabel('sigma amplitude');
title('sigma amplitude');
grid on;

%% comparison of stress baseline and amplitude for TGF beta

TGF_values = [1, 1];
TGFF_values = [1, 5];
colors = ['b', 'r']; 
labels = {'No TGF', 'TGF'};

baseline_results = zeros(length(TGFF_values), length(epsilon_ext_values));
amplitude_results = zeros(length(TGFF_values), length(epsilon_ext_values));

for idx = 1:length(TGFF_values)
    TGF_current = TGF_values(idx);
    TGFF_current = TGFF_values(idx);
    E_t0f_current = TGFF_values(idx) * E_t0f;
    sigma_min_vals = zeros(size(epsilon_ext_values));
    sigma_max_vals = zeros(size(epsilon_ext_values));

    for j = 1:length(epsilon_ext_values)
        epsilon_ext_current = epsilon_ext_values(j);
        sigma_vals = [];
        t_vals = [];

        %  initial conditions
        epsilon_yita_prev = 0;
        initial_guess = zeros(1, 15);
        initial_guess(8) = epsilon_yita_prev;
        initial_guess(10) = C0_base;

        for t = time_range
            options = optimoptions('fsolve', 'Display', 'none');
            sol = fsolve(@(variables) system_of_equations(variables, t, epsilon_yita_prev, ...
                C0_base, amplitude, E_p, E_t0, E_c, beta, alpha, yita, n, ...
                epsilon_ext_current, dt, l, m, TGF_current, ...
                E_t0f_current, E_cf, betaf, alphaf, Rho0, TGFF_current, lf, mf, ...
                E_col0, lcol, mcol), initial_guess, options);

            sigma_vals(end + 1) = sol(5);
            t_vals(end + 1) = t;
            epsilon_yita_prev = sol(8);
            initial_guess = sol;
        end

        % extract only steady-state values
        steady_inds = t_vals >= time_steady;
        sigma_min_vals(j) = min(sigma_vals(steady_inds));
        sigma_max_vals(j) = max(sigma_vals(steady_inds));
    end

    baseline_results(idx, :) = sigma_min_vals;
    amplitude_results(idx, :) = sigma_max_vals - sigma_min_vals;
end

%  baseline stress 
figure;
hold on;
for idx = 1:length(TGFF_values)
    plot(epsilon_ext_values, baseline_results(idx, :), [colors(idx) '-o'], 'LineWidth', 2, 'DisplayName', labels{idx});
end
xlabel('Epsilon_{ext}');
ylabel('Baseline sigma');
title('Baseline sigma');
legend;
grid on;

%  stress amplitude 
figure;
hold on;
for idx = 1:length(TGFF_values)
    plot(epsilon_ext_values, amplitude_results(idx, :), [colors(idx) '-o'], 'LineWidth', 2, 'DisplayName', labels{idx});
end
xlabel('Epsilon_{ext}');
ylabel('stress amplitude');
title('stress amplitude for TGF beta');
legend;
grid on;


%% comparison of contractility baseline and amplitude for TGF beta

TGF_values = [1, 1];
TGFF_values = [1, 5];
colors = ['b', 'r'];
labels = {'No TGF', 'TGF'};

baseline_C_results = zeros(length(TGFF_values), length(epsilon_ext_values));
amplitude_C_results = zeros(length(TGFF_values), length(epsilon_ext_values));

for idx = 1:length(TGFF_values)
    TGF_current = TGF_values(idx);
    TGFF_current = TGFF_values(idx);
    E_t0f_current = TGFF_current * E_t0f;
    C_min_vals = zeros(size(epsilon_ext_values));
    C_max_vals = zeros(size(epsilon_ext_values));

    for j = 1:length(epsilon_ext_values)
        epsilon_ext_current = epsilon_ext_values(j);
        C_vals = [];
        t_vals = [];

        %  initial conditions
        epsilon_yita_prev = 0;
        initial_guess = zeros(1, 15);
        initial_guess(8) = epsilon_yita_prev;
        initial_guess(10) = C0_base;

        for t = time_range
            options = optimoptions('fsolve', 'Display', 'none');
            sol = fsolve(@(variables) system_of_equations(variables, t, epsilon_yita_prev, ...
                C0_base, amplitude, E_p, E_t0, E_c, beta, alpha, yita, n, ...
                epsilon_ext_current, dt, l, m, TGF_current, ...
                E_t0f_current, E_cf, betaf, alphaf, Rho0, TGFF_current, lf, mf, ...
                E_col0, lcol, mcol), initial_guess, options);

            C_vals(end + 1) = sol(4);  % contractility C
            t_vals(end + 1) = t;
            epsilon_yita_prev = sol(8);
            initial_guess = sol;
        end

        % extract steady-state values
        steady_inds = t_vals >= time_steady;
        C_min_vals(j) = min(C_vals(steady_inds));
        C_max_vals(j) = max(C_vals(steady_inds));
    end

    baseline_C_results(idx, :) = C_min_vals;
    amplitude_C_results(idx, :) = C_max_vals - C_min_vals;
end

% baseline contractility
figure;
hold on;
for idx = 1:length(TGFF_values)
    plot(epsilon_ext_values, baseline_C_results(idx, :), [colors(idx) '-o'], 'LineWidth', 2, 'DisplayName', labels{idx});
end
xlabel('Epsilon_{ext}');
ylabel('Baseline contractility');
title('Baseline contractility for TGF beta');
legend;
grid on;

% contractility amplitude
figure;
hold on;
for idx = 1:length(TGFF_values)
    plot(epsilon_ext_values, amplitude_C_results(idx, :), [colors(idx) '-o'], 'LineWidth', 2, 'DisplayName', labels{idx});
end
xlabel('Epsilon_{ext}');
ylabel('Contractility amplitude');
title('Contractility amplitude for TGF beta');
legend;
grid on;



%% stress baseline and amplitude for TGF beta without strain stiffening
l_nostiff = 0;
TGF_values = [1, 1];
TGFF_values = [1, 5];
colors = ['b', 'r']; 
labels = {'No TGF', 'TGF'};

baseline_results = zeros(length(TGFF_values), length(epsilon_ext_values));
amplitude_results = zeros(length(TGFF_values), length(epsilon_ext_values));

for idx = 1:length(TGFF_values)
    TGF_current = TGF_values(idx);
    TGFF_current = TGFF_values(idx);
    E_t0f_current = TGFF_values(idx) * E_t0f;
    sigma_min_vals = zeros(size(epsilon_ext_values));
    sigma_max_vals = zeros(size(epsilon_ext_values));

    for j = 1:length(epsilon_ext_values)
        epsilon_ext_current = epsilon_ext_values(j);
        sigma_vals = [];
        t_vals = [];

        %  initial conditions
        epsilon_yita_prev = 0;
        initial_guess = zeros(1, 15);
        initial_guess(8) = epsilon_yita_prev;
        initial_guess(10) = C0_base;

        for t = time_range
            options = optimoptions('fsolve', 'Display', 'none');
            sol = fsolve(@(variables) system_of_equations(variables, t, epsilon_yita_prev, ...
                C0_base, amplitude, E_p, E_t0, E_c, beta, alpha, yita, n, ...
                epsilon_ext_current, dt, l_nostiff, m, TGF_current, ...
                E_t0f_current, E_cf, betaf, alphaf, Rho0, TGFF_current, lf, mf, ...
                E_col0, lcol, mcol), initial_guess, options);

            sigma_vals(end + 1) = sol(5);
            t_vals(end + 1) = t;
            epsilon_yita_prev = sol(8);
            initial_guess = sol;
        end

        % extract only steady-state values
        steady_inds = t_vals >= time_steady;
        sigma_min_vals(j) = min(sigma_vals(steady_inds));
        sigma_max_vals(j) = max(sigma_vals(steady_inds));
    end

    baseline_results(idx, :) = sigma_min_vals;
    amplitude_results(idx, :) = sigma_max_vals - sigma_min_vals;
end

%  baseline force 
figure;
hold on;
for idx = 1:length(TGFF_values)
    plot(epsilon_ext_values, baseline_results(idx, :), [colors(idx) '-o'], 'LineWidth', 2, 'DisplayName', labels{idx});
end
xlabel('Epsilon_{ext}');
ylabel('Baseline stress');
title('Baseline stress for TGF beta without strain stiffening');
legend;
grid on;

%  cardiac force 
figure;
hold on;
for idx = 1:length(TGFF_values)
    plot(epsilon_ext_values, amplitude_results(idx, :), [colors(idx) '-o'], 'LineWidth', 2, 'DisplayName', labels{idx});
end
xlabel('Epsilon_{ext}');
ylabel('stress amplitude');
title('stress amplitude for TGF beta without strain stiffening');
legend;
grid on;


%% contractility baseline and amplitude for TGF beta without strain stiffening
l_nostiff = 0;
TGF_values = [1, 1];
TGFF_values = [1, 5];
colors = ['b', 'r'];
labels = {'No TGF', 'TGF'};

baseline_C_results = zeros(length(TGFF_values), length(epsilon_ext_values));
amplitude_C_results = zeros(length(TGFF_values), length(epsilon_ext_values));

for idx = 1:length(TGFF_values)
    TGF_current = TGF_values(idx);
    TGFF_current = TGFF_values(idx);
    E_t0f_current = TGFF_current * E_t0f;
    C_min_vals = zeros(size(epsilon_ext_values));
    C_max_vals = zeros(size(epsilon_ext_values));

    for j = 1:length(epsilon_ext_values)
        epsilon_ext_current = epsilon_ext_values(j);
        C_vals = [];
        t_vals = [];

        %  initial conditions
        epsilon_yita_prev = 0;
        initial_guess = zeros(1, 15);
        initial_guess(8) = epsilon_yita_prev;
        initial_guess(10) = C0_base;

        for t = time_range
            options = optimoptions('fsolve', 'Display', 'none');
            sol = fsolve(@(variables) system_of_equations(variables, t, epsilon_yita_prev, ...
                C0_base, amplitude, E_p, E_t0, E_c, beta, alpha, yita, n, ...
                epsilon_ext_current, dt, l_nostiff, m, TGF_current, ...
                E_t0f_current, E_cf, betaf, alphaf, Rho0, TGFF_current, lf, mf, ...
                E_col0, lcol, mcol), initial_guess, options);

            C_vals(end + 1) = sol(4);  % contractility
            t_vals(end + 1) = t;
            epsilon_yita_prev = sol(8);
            initial_guess = sol;
        end

        % extract only steady-state values
        steady_inds = t_vals >= time_steady;
        C_min_vals(j) = min(C_vals(steady_inds));
        C_max_vals(j) = max(C_vals(steady_inds));
    end

    baseline_C_results(idx, :) = C_min_vals;
    amplitude_C_results(idx, :) = C_max_vals - C_min_vals;
end

% baseline contractility 
figure;
hold on;
for idx = 1:length(TGFF_values)
    plot(epsilon_ext_values, baseline_C_results(idx, :), [colors(idx) '-o'], ...
        'LineWidth', 2, 'DisplayName', labels{idx});
end
xlabel('Epsilon_{ext}');
ylabel('Baseline contractility');
title('Baseline contractility for TGF beta (no strain stiffening)');
legend;
grid on;

% cardiac contractility amplitude 
figure;
hold on;
for idx = 1:length(TGFF_values)
    plot(epsilon_ext_values, amplitude_C_results(idx, :), [colors(idx) '-o'], ...
        'LineWidth', 2, 'DisplayName', labels{idx});
end
xlabel('Epsilon_{ext}');
ylabel('Contractility amplitude');
title('Contractility amplitude for TGF beta (no strain stiffening)');
legend;
grid on;



%% stress baseline and amplitude for IVS

TGFF_current = 5;
E_t0f_treated = TGFF_current * E_t0f;
IVS_dose = [1, 0.15, 0.1, 0.05];
colors = lines(length(IVS_dose)); 
labels = {'IVS 0', 'IVS 1', 'IVS 3', 'IVS 10'};

baseline_results = zeros(length(IVS_dose), length(epsilon_ext_values));
amplitude_results = zeros(length(IVS_dose), length(epsilon_ext_values));

for idx = 1:length(IVS_dose)
    E_t0f_current = IVS_dose(idx) * E_t0f_treated;
    sigma_min_vals = zeros(size(epsilon_ext_values));
    sigma_max_vals = zeros(size(epsilon_ext_values));

    for j = 1:length(epsilon_ext_values)
        epsilon_ext_current = epsilon_ext_values(j);
        sigma_vals = [];
        t_vals = [];

        %  initial conditions
        epsilon_yita_prev = 0;
        initial_guess = zeros(1, 15);
        initial_guess(8) = epsilon_yita_prev;
        initial_guess(10) = C0_base;

        for t = time_range
            options = optimoptions('fsolve', 'Display', 'none');
            sol = fsolve(@(variables) system_of_equations(variables, t, epsilon_yita_prev, ...
                C0_base, amplitude, E_p, E_t0, E_c, beta, alpha, yita, n, ...
                epsilon_ext_current, dt, l, m, TGF, ...
                E_t0f_current, E_cf, betaf, alphaf, Rho0, TGFF_current, lf, mf, ...
                E_col0, lcol, mcol), initial_guess, options);

            sigma_vals(end + 1) = sol(5);
            t_vals(end + 1) = t;
            epsilon_yita_prev = sol(8);
            initial_guess = sol;
        end

        % extract only steady-state values
        steady_inds = t_vals >= time_steady;
        sigma_min_vals(j) = min(sigma_vals(steady_inds));
        sigma_max_vals(j) = max(sigma_vals(steady_inds));
    end

    baseline_results(idx, :) = sigma_min_vals;
    amplitude_results(idx, :) = sigma_max_vals - sigma_min_vals;
end


%  baseline force 
figure;
hold on;
for idx = 1:length(IVS_dose)
    plot(epsilon_ext_values, baseline_results(idx, :), '-o', ...
    'Color', colors(idx, :), 'LineWidth', 2, 'MarkerSize', 6, ...
    'DisplayName', labels{idx});
end
xlabel('Epsilon_{ext}');
ylabel('Baseline sigma');
title('Baseline sigma for IVS');
legend;
grid on;

%  sigma amplitude
figure;
hold on;
for idx = 1:length(IVS_dose)
    plot(epsilon_ext_values, amplitude_results(idx, :), '-o', ...
    'Color', colors(idx, :), 'LineWidth', 2, 'MarkerSize', 6, ...
    'DisplayName', labels{idx});
end
xlabel('Epsilon_{ext}');
ylabel('sigma amplitude');
title('sigma amplitude for IVS');
legend;
grid on;


%% contractility baseline and cardiac contractility for IVS

TGFF_current = 5;
E_t0f_treated = TGFF_current * E_t0f;
IVS_dose = [1, 0.15, 0.1, 0.05];
colors = lines(length(IVS_dose));
labels = {'IVS 0', 'IVS 1', 'IVS 3', 'IVS 10'};

baseline_results_C = zeros(length(IVS_dose), length(epsilon_ext_values));
amplitude_results_C = zeros(length(IVS_dose), length(epsilon_ext_values));

for idx = 1:length(IVS_dose)
    E_t0f_current = IVS_dose(idx) * E_t0f_treated;
    C_min_vals = zeros(size(epsilon_ext_values));
    C_max_vals = zeros(size(epsilon_ext_values));

    for j = 1:length(epsilon_ext_values)
        epsilon_ext_current = epsilon_ext_values(j);
        C_vals = [];
        t_vals = [];

        epsilon_yita_prev = 0;
        initial_guess = zeros(1, 15);
        initial_guess(8) = epsilon_yita_prev;
        initial_guess(10) = C0_base;

        for t = time_range
            options = optimoptions('fsolve', 'Display', 'none');
            sol = fsolve(@(variables) system_of_equations(variables, t, epsilon_yita_prev, ...
                C0_base, amplitude, E_p, E_t0, E_c, beta, alpha, yita, n, ...
                epsilon_ext_current, dt, l, m, TGF, ...
                E_t0f_current, E_cf, betaf, alphaf, Rho0, TGFF_current, lf, mf, ...
                E_col0, lcol, mcol), initial_guess, options);

            C_vals(end + 1) = sol(4);  % contractility
            t_vals(end + 1) = t;
            epsilon_yita_prev = sol(8);
            initial_guess = sol;
        end

        % steady-state 
        steady_inds = t_vals >= time_steady;
        C_min_vals(j) = min(C_vals(steady_inds));
        C_max_vals(j) = max(C_vals(steady_inds));
    end

    baseline_results_C(idx, :) = C_min_vals;
    amplitude_results_C(idx, :) = C_max_vals - C_min_vals;
end

%  baseline contractility 
figure;
hold on;
for idx = 1:length(IVS_dose)
    plot(epsilon_ext_values, baseline_results_C(idx, :), '-o', ...
        'Color', colors(idx, :), 'LineWidth', 2, 'MarkerSize', 6, ...
        'DisplayName', labels{idx});
end
xlabel('Epsilon_{ext}');
ylabel('Baseline contractility');
title('Baseline contractility for IVS');
legend;
grid on;

% contractility amplitude
figure;
hold on;
for idx = 1:length(IVS_dose)
    plot(epsilon_ext_values, amplitude_results_C(idx, :), '-o', ...
        'Color', colors(idx, :), 'LineWidth', 2, 'MarkerSize', 6, ...
        'DisplayName', labels{idx});
end
xlabel('Epsilon_{ext}');
ylabel('contractility amplitude');
title('contractility amplitude for IVS');
legend;
grid on;


%%  contractility and stress vs. time (no stretching)

% time settings
time_plot = linspace(0, 20, 10000);
dt_plot = time_plot(2) - time_plot(1);

epsilon_ext_current = 0; % constant stretch
epsilon_yita_prev = 0;
initial_guess = zeros(1, 15);
initial_guess(8) = epsilon_yita_prev;
initial_guess(10) = C0_base;

C_values = [];
Sigma_values = [];
time_log = [];

for t = time_plot
    options = optimoptions('fsolve', 'Display', 'none');
    sol = fsolve(@(variables) system_of_equations(variables, t, epsilon_yita_prev, ...
        C0_base, amplitude, E_p, E_t0, E_c, beta, alpha, yita, n, ...
        epsilon_ext_current, dt_plot, l, m, TGF, ...
        E_t0f, E_cf, betaf, alphaf, Rho0, TGFF, lf, mf, ...
        E_col0, lcol, mcol), initial_guess, options);

    C_values(end + 1) = sol(4);      % C
    Sigma_values(end + 1) = sol(5);  % sigma
    time_log(end + 1) = t;

    epsilon_yita_prev = sol(8);
    initial_guess = sol;
end

% only plot between 10 and 20 sec
mask = time_log >= 10;
t_plot = time_log(mask);
C_plot = C_values(mask);
Sigma_plot = Sigma_values(mask);

% plot C 
figure;
plot(t_plot, C_plot, 'b-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Contractility (C)');
title('Contractility vs. Time');
grid on;

% plot sigma
figure;
plot(t_plot, Sigma_plot, 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Stress');
title('Stress vs. Time');
grid on;


%% Stress vs. Time (with and without active contractility)

Sigma_with_C = Sigma_values;

C0_base_passive = 0;
amplitude_passive = 0;

Sigma_no_C = [];
epsilon_yita_prev = 0;
initial_guess = zeros(1, 15);
initial_guess(8) = epsilon_yita_prev;
initial_guess(10) = C0_base_passive;

for t = time_plot
    options = optimoptions('fsolve', 'Display', 'none');
    sol = fsolve(@(variables) system_of_equations(variables, t, epsilon_yita_prev, ...
        C0_base_passive, amplitude_passive, E_p, E_t0, E_c, beta, alpha, yita, n, ...
        epsilon_ext_current, dt_plot, l, m, TGF, ...
        E_t0f, E_cf, betaf, alphaf, Rho0, TGFF, lf, mf, ...
        E_col0, lcol, mcol), initial_guess, options);

    Sigma_no_C(end + 1) = sol(5);
    epsilon_yita_prev = sol(8);
    initial_guess = sol;
end

% only plot for 10 to 20 sec
Sigma_with_plot = Sigma_with_C(mask);
Sigma_no_plot = Sigma_no_C(mask);

figure;
plot(t_plot, Sigma_with_plot, 'r-', 'LineWidth', 2); hold on;
plot(t_plot, Sigma_no_plot, 'k--', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Stress');
title('Stress vs. time with and without active CM contractility');
legend('with contractility', 'without contractility (C₀ = 0, amp = 0)');
grid on;


%% contractility and stress vs. Time for TGF beta

epsilon_ext_current = 0; % constant stretch
TGFF_current = 5;
TGFF_values = [TGFF, TGFF_current];
colors = {'b', 'r'}; 
labels = {'No TGF', 'TGF'};

figure_C = figure; hold on;
figure_Sigma = figure; hold on;

for i = 1:length(TGFF_values)
    TGFF_current = TGFF_values(i);
    E_t0f_current = TGFF_values(i) * E_t0f;
    epsilon_yita_prev = 0;
    initial_guess = zeros(1, 15);
    initial_guess(8) = epsilon_yita_prev;
    initial_guess(10) = C0_base;

    C_values = [];
    Sigma_values = [];
    time_log = [];

    for t = time_plot
        options = optimoptions('fsolve', 'Display', 'none');
        sol = fsolve(@(variables) system_of_equations(variables, t, epsilon_yita_prev, ...
            C0_base, amplitude, E_p, E_t0, E_c, beta, alpha, yita, n, ...
            epsilon_ext_current, dt_plot, l, m, TGF, ...
            E_t0f_current, E_cf, betaf, alphaf, Rho0, TGFF_current, lf, mf, ...
            E_col0, lcol, mcol), initial_guess, options);

        C_values(end + 1) = sol(4);      % C
        Sigma_values(end + 1) = sol(5);  % sigma
        time_log(end + 1) = t;

        epsilon_yita_prev = sol(8);
        initial_guess = sol;
    end

    % only plot values between 10 and 20 sec
    mask = time_log >= 10;
    t_plot = time_log(mask);
    C_plot = C_values(mask);
    Sigma_plot = Sigma_values(mask);

    figure(figure_C);
    plot(t_plot, C_plot, 'Color', colors{i}, 'LineWidth', 2, 'DisplayName', labels{i});

    figure(figure_Sigma);
    plot(t_plot, Sigma_plot, 'Color', colors{i}, 'LineWidth', 2, 'DisplayName', labels{i});
end

% contractility
figure(figure_C);
xlabel('Time (s)');
ylabel('Contractility (C)');
title('Contractility vs. Time  for TGF beta');
legend;
grid on;

% stress
figure(figure_Sigma);
xlabel('Time (s)');
ylabel('Stress');
title('Stress vs. Time for TGF beta');
legend;
grid on;



%% Contractility and stress vs. Time for IVS

epsilon_ext_current = 0; % constant stretch
TGFF_current = 5;
E_t0f_treated = TGFF_current * E_t0f;
IVS_dose = [1, 0.15, 0.1, 0.05];
colors = {'k', 'r', 'g', 'b'}; 
labels = {'IVS 0', 'IVS 1', 'IVS 3', 'IVS 10'};

figure_CIVS = figure; hold on;
figure_SigmaIVS = figure; hold on;

for i = 1:length(IVS_dose)
    E_t0f_current = IVS_dose(i) * E_t0f_treated;
    epsilon_yita_prev = 0;
    initial_guess = zeros(1, 15);
    initial_guess(8) = epsilon_yita_prev;
    initial_guess(10) = C0_base;

    C_values = [];
    Sigma_values = [];
    time_log = [];

    for t = time_plot
        options = optimoptions('fsolve', 'Display', 'none');
        sol = fsolve(@(variables) system_of_equations(variables, t, epsilon_yita_prev, ...
            C0_base, amplitude, E_p, E_t0, E_c, beta, alpha, yita, n, ...
            epsilon_ext_current, dt_plot, l, m, TGF, ...
            E_t0f_current, E_cf, betaf, alphaf, Rho0, TGFF_current, lf, mf, ...
            E_col0, lcol, mcol), initial_guess, options);

        C_values(end + 1) = sol(4);      % C
        Sigma_values(end + 1) = sol(5);  % sigma
        time_log(end + 1) = t;

        epsilon_yita_prev = sol(8);
        initial_guess = sol;
    end

    % only values between 10 and 20 sec
    mask = time_log >= 10;
    t_plot = time_log(mask);
    C_plot = C_values(mask);
    Sigma_plot = Sigma_values(mask);

    figure(figure_CIVS);
    plot(t_plot, C_plot, 'Color', colors{i}, 'LineWidth', 2, 'DisplayName', labels{i});

    figure(figure_SigmaIVS);
    plot(t_plot, Sigma_plot, 'Color', colors{i}, 'LineWidth', 2, 'DisplayName', labels{i});
end

% contractility
figure(figure_CIVS);
xlabel('Time (s)');
ylabel('Contractility (C)');
title('Contractility vs. Time for IVS');
legend;
grid on;

% stress
figure(figure_SigmaIVS);
xlabel('Time (s)');
ylabel('Stress');
title('Stress vs. Time for IVS');
legend;
grid on;

%%  Contractility and stress vs. post stiffness (Ep)

Ep_values = 0:1:50;  % from 0 to 50 in unit steps
C_min = zeros(size(Ep_values));
C_avg = zeros(size(Ep_values));
C_max = zeros(size(Ep_values));
Sigma_min = zeros(size(Ep_values));
Sigma_avg = zeros(size(Ep_values));
Sigma_max = zeros(size(Ep_values));

time_range = time_range; 
dt_plot = time_range(2) - time_range(1);
epsilon_ext_current = 0;  % fixed stretch

for k = 1:length(Ep_values)
    Ep = Ep_values(k);

    C_values = [];
    Sigma_values = [];
    time_log = [];

    epsilon_yita_prev = 0;
    initial_guess = zeros(1, 15);
    initial_guess(8) = epsilon_yita_prev;
    initial_guess(10) = C0_base;

    for t = time_range
        options = optimoptions('fsolve', 'Display', 'none');
        sol = fsolve(@(variables) system_of_equations(variables, t, epsilon_yita_prev, ...
            C0_base, amplitude, Ep, E_t0, E_c, beta, alpha, yita, n, ...
            epsilon_ext_current, dt_plot, l, m, TGF, ...
            E_t0f, E_cf, betaf, alphaf, Rho0, TGFF, lf, mf, ...
            E_col0, lcol, mcol), initial_guess, options);

        C_values(end + 1) = sol(4);       % contractility
        Sigma_values(end + 1) = sol(5);   % stress
        time_log(end + 1) = t;

        epsilon_yita_prev = sol(8);
        initial_guess = sol;
    end

    % extract values from steady
    mask = time_log >= time_steady;
    C_steady = C_values(mask);
    Sigma_steady = Sigma_values(mask);

    C_min(k) = min(C_steady);
    C_avg(k) = mean(C_steady);
    C_max(k) = max(C_steady);
    Sigma_min(k) = min(Sigma_steady);
    Sigma_avg(k) = mean(Sigma_steady);
    Sigma_max(k) = max(Sigma_steady);
end

% contractility vs Ep
figure;
plot(Ep_values, C_min, 'b--', Ep_values, C_avg, 'b-', Ep_values, C_max, 'b-.', 'LineWidth', 2);
xlabel('Post Stiffness (E_p)');
ylabel('Contractility (C)');
title('Contractility vs. Post Stiffness (E_p)');
legend('Min (baseline)', 'Mean', 'Max (peak)', 'Location', 'best');
grid on;

% stress vs Ep
figure;
plot(Ep_values, Sigma_min, 'r--', Ep_values, Sigma_avg, 'r-', Ep_values, Sigma_max, 'r-.', 'LineWidth', 2);
xlabel('Post Stiffness (E_p)');
ylabel('Stress');
title('Stress vs. Post Stiffness (E_p)');
legend('Min (baseline)', 'Mean', 'Max (peak)', 'Location', 'best');
grid on;

%%  base and amplitude of contractility and stress vs. post stiffness (Ep) with and without strain stiffening

Ep_values = 0:1:50; 
C_min_noss = zeros(size(Ep_values));
C_avg_noss = zeros(size(Ep_values));
C_max_noss = zeros(size(Ep_values));
Sigma_min_noss = zeros(size(Ep_values));
Sigma_avg_noss = zeros(size(Ep_values));
Sigma_max_noss = zeros(size(Ep_values));

time_range = time_range;  
dt_plot = time_range(2) - time_range(1);
epsilon_ext_current = 0;  % fixed stretch
l_nostiff = 0;

for k = 1:length(Ep_values)
    Ep = Ep_values(k);

    C_values_noss = [];
    Sigma_values_noss = [];
    time_log = [];

    epsilon_yita_prev = 0;
    initial_guess = zeros(1, 15);
    initial_guess(8) = epsilon_yita_prev;
    initial_guess(10) = C0_base;

    for t = time_range
        options = optimoptions('fsolve', 'Display', 'none');
        sol = fsolve(@(variables) system_of_equations(variables, t, epsilon_yita_prev, ...
            C0_base, amplitude, Ep, E_t0, E_c, beta, alpha, yita, n, ...
            epsilon_ext_current, dt_plot, l_nostiff, m, TGF, ...
            E_t0f, E_cf, betaf, alphaf, Rho0, TGFF, lf, mf, ...
            E_col0, lcol, mcol), initial_guess, options);

        C_values_noss(end + 1) = sol(4);       % contractility
        Sigma_values_noss(end + 1) = sol(5);   % stress
        time_log(end + 1) = t;

        epsilon_yita_prev = sol(8);
        initial_guess = sol;
    end

    % extract values from steady
    mask = time_log >= time_steady;
    C_steady_noss = C_values_noss(mask);
    Sigma_steady_noss = Sigma_values_noss(mask);

    C_min_noss(k) = min(C_steady_noss);
    C_avg_noss(k) = mean(C_steady_noss);
    C_max_noss(k) = max(C_steady_noss);
    Sigma_min_noss(k) = min(Sigma_steady_noss);
    Sigma_avg_noss(k) = mean(Sigma_steady_noss);
    Sigma_max_noss(k) = max(Sigma_steady_noss);
end

% base contractility
figure;
plot(Ep_values, C_min, 'b-', 'LineWidth', 2); hold on;
plot(Ep_values, C_min_noss, 'b--', 'LineWidth', 2);
xlabel('Post stiffness (E_p)');
ylabel('Baseline contractility');
title('Baseline contractility vs. E_p (with and without stiffening)');
legend('With stiffening', 'Without stiffening');
grid on;

% contractility amplitude
C_amp = C_max - C_min;
C_amp_noss = C_max_noss - C_min_noss;

figure;
plot(Ep_values, C_amp, 'b-', 'LineWidth', 2); hold on;
plot(Ep_values, C_amp_noss, 'b--', 'LineWidth', 2);
xlabel('Post stiffness (E_p)');
ylabel('Contractility amplitude');
title('Contractility amplitude vs. E_p with and without stiffening');
legend('With stiffening', 'Without stiffening');
grid on;

% stress baseline
figure;
plot(Ep_values, Sigma_min, 'r-', 'LineWidth', 2); hold on;
plot(Ep_values, Sigma_min_noss, 'r--', 'LineWidth', 2);
xlabel('Post stiffness (E_p)');
ylabel('Baseline stress');
title('Baseline stress vs. E_p with and without stiffening');
legend('With stiffening', 'Without stiffening');
grid on;

% stress amplitude
Sigma_amp = Sigma_max - Sigma_min;
Sigma_amp_noss = Sigma_max_noss - Sigma_min_noss;

figure;
plot(Ep_values, Sigma_amp, 'r-', 'LineWidth', 2); hold on;
plot(Ep_values, Sigma_amp_noss, 'r--', 'LineWidth', 2);
xlabel('Post stiffness (E_p)');
ylabel('Stress amplitude');
title('Stress amplitude vs. E_p with and without stiffening');
legend('With stiffening', 'Without stiffening');
grid on;




%% Stress vs. Et strain

epsilont_range = linspace(0, 0.3, 100);  
Sigmat = E_t0 * epsilont_range + (E_t0 * l) ./ (m + 1) .* epsilont_range.^(m + 1);  % F(6) formulation

figure;
plot(epsilont_range, Sigmat, 'k-', 'LineWidth', 2);
xlabel('strain of cyto in tension (\epsilon_t)');
ylabel('stress with stiffening');
title('Stress vs. \epsilon_t');
grid on;


% no stiffening parameters
Sigmat_noss = E_t0 * epsilont_range;

figure;
plot(epsilont_range, Sigmat_noss, 'k-', 'LineWidth', 2);
xlabel('strain (\epsilon_t)');
ylabel('stress without stiffening');
title('Stress vs. \epsilon_t without strain stiffening');
grid on;

%%  Contractility and stress vs. time with and without stiffening

C_with_stiff = [];
Sigma_with_stiff = [];

C_no_stiff = [];
Sigma_no_stiff = [];

epsilon_ext_current = 0; % constant stretch
time_log = time_plot;

% with strain stiffening
epsilon_yita_prev = 0;
initial_guess = zeros(1, 15);
initial_guess(8) = epsilon_yita_prev;
initial_guess(10) = C0_base;

for t = time_plot
    options = optimoptions('fsolve', 'Display', 'none');
    sol = fsolve(@(variables) system_of_equations(variables, t, epsilon_yita_prev, ...
        C0_base, amplitude, E_p, E_t0, E_c, beta, alpha, yita, n, ...
        epsilon_ext_current, dt_plot, l, m, TGF, ...
        E_t0f, E_cf, betaf, alphaf, Rho0, TGFF, lf, mf, ...
        E_col0, lcol, mcol), initial_guess, options);

    C_with_stiff(end + 1) = sol(4);
    Sigma_with_stiff(end + 1) = sol(5);
    epsilon_yita_prev = sol(8);
    initial_guess = sol;
end

% Without strain stiffening
epsilon_yita_prev = 0;
initial_guess = zeros(1, 15);
initial_guess(8) = epsilon_yita_prev;
initial_guess(10) = C0_base;
l_nostiff = 0;

for t = time_plot
    options = optimoptions('fsolve', 'Display', 'none');
    sol = fsolve(@(variables) system_of_equations(variables, t, epsilon_yita_prev, ...
        C0_base, amplitude, E_p, E_t0, E_c, beta, alpha, yita, n, ...
        epsilon_ext_current, dt_plot, l_nostiff, m, TGF, ...
        E_t0f, E_cf, betaf, alphaf, Rho0, TGFF, lf, mf, ...
        E_col0, lcol, mcol), initial_guess, options);

    C_no_stiff(end + 1) = sol(4);
    Sigma_no_stiff(end + 1) = sol(5);
    epsilon_yita_prev = sol(8);
    initial_guess = sol;
end

% Only plot from 10 to 20 sec
mask = time_log >= 10;
t_plot = time_log(mask);

C_with = C_with_stiff(mask);
C_no = C_no_stiff(mask);
Sigma_with = Sigma_with_stiff(mask);
Sigma_no = Sigma_no_stiff(mask);

% C
figure;
plot(t_plot, C_with, 'b-', 'LineWidth', 2); hold on;
plot(t_plot, C_no, 'b--', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Contractility (C)');
title('Contractility vs. Time with and without stiffening');
legend('With strain stiffening', 'without strain stiffening');
grid on;

% sigma
figure;
plot(t_plot, Sigma_with, 'r-', 'LineWidth', 2); hold on;
plot(t_plot, Sigma_no, 'r--', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Stress');
title('Stress vs. Time with and without stiffening');
legend('With strain stiffening', 'without strain stiffening');
grid on;



%% Pillar strain vs. post stiffness (E_p)

Ep_values = 0:1:50; 
epsilon_p_min = zeros(size(Ep_values));
epsilon_p_max = zeros(size(Ep_values));
epsilon_p_avg = zeros(size(Ep_values));

time_plot = time_range;
dt_plot = time_plot(2) - time_plot(1);
epsilon_ext_current = 0;  % fixed stretch

for k = 1:length(Ep_values)
    Ep = Ep_values(k);

    epsilon_p_vals = [];
    time_log = [];
    epsilon_yita_prev = 0;
    initial_guess = zeros(1, 15);
    initial_guess(8) = epsilon_yita_prev;
    initial_guess(10) = C0_base;

    for t = time_plot
        options = optimoptions('fsolve', 'Display', 'none');
        sol = fsolve(@(variables) system_of_equations(variables, t, epsilon_yita_prev, ...
            C0_base, amplitude, Ep, E_t0, E_c, beta, alpha, yita, n, ...
            epsilon_ext_current, dt_plot, l, m, TGF, ...
            E_t0f, E_cf, betaf, alphaf, Rho0, TGFF, lf, mf, ...
            E_col0, lcol, mcol), initial_guess, options);

        epsilon_p_vals(end + 1) = sol(1);  % ε_p
        time_log(end + 1) = t;
        epsilon_yita_prev = sol(8);
        initial_guess = sol;
    end

    % Only time 10 to 20 s
    mask = time_log >= time_steady;
    eps_p_steady = epsilon_p_vals(mask);
    epsilon_p_min(k) = min(eps_p_steady);
    epsilon_p_max(k) = max(eps_p_steady);
    epsilon_p_avg(k) = mean(eps_p_steady);
end

figure;
plot(Ep_values, epsilon_p_min, 'b--', 'LineWidth', 2); hold on;
plot(Ep_values, epsilon_p_avg, 'b-', 'LineWidth', 2);
plot(Ep_values, epsilon_p_max, 'b-.', 'LineWidth', 2);
xlabel('Post stiffness (E_p)');
ylabel('Pillar strain (\epsilon_p)');
title('Pillar strain vs. post stiffness (E_p)');
legend('Min', 'Mean', 'Max', 'Location', 'northeast');
grid on;


relative_shortening = abs(epsilon_p_max - epsilon_p_min);

figure;
plot(Ep_values, relative_shortening, 'k-', 'LineWidth', 2);
xlabel('Post stiffness (E_p)');
ylabel('cell shortening');
title('cell shortening vs E_p');
grid on;


%% function
function F = system_of_equations(variables, time_step, epsilon_yita_prev, C0_base, amplitude, E_p, E_t0, E_c, beta, alpha, yita, n, epsilon_ext, dt, l, m, TGF, E_t0f, E_cf, betaf, alphaf, Rho0, TGFF, lf, mf, E_col0, lcol, mcol)
    epsilon_p = variables(1);
    epsilon_t = variables(2);
    epsilon_c = variables(3);
    C = variables(4);
    Sigma = variables(5);
    Sigma_yita = variables(6);
    Sigma_t = variables(7);
    epsilon_yita = variables(8);
    C0 = variables(9);
    epsilon_cf = variables(10);
    Sigma_tf = variables(11);
    Rho = variables(12);
    epsilon_tf = variables(13);
    epsilon_col = variables(14);
    Sigma_col = variables(15);

    
    % update C0
    C0_current = TGF * (C0_base + amplitude * ((abs(sin(pi * time_step))) .^ 5));

    Rho0 = TGFF * Rho0;

    % system of equations
    F(1) = epsilon_c + beta * (C - C0_current) - alpha * ((Sigma_t) .^ n);
    F(2) = Sigma - E_p * epsilon_p;
    F(3) = Sigma - Sigma_t;
    F(4) = Sigma - (Sigma_col + Sigma_yita);
    F(5) = (epsilon_yita - epsilon_yita_prev) - (Sigma_yita / yita) * dt;
    F(6) = Sigma_t - E_t0 * epsilon_t - ((E_t0 * l) / (m + 1)) * epsilon_t ^ (m + 1);
    F(7) = Sigma_t - (C + E_c * epsilon_c);
    F(8) = epsilon_yita - epsilon_col;
    F(9) = 2 * epsilon_p + epsilon_c + epsilon_t + epsilon_cf + epsilon_tf + epsilon_yita - epsilon_ext;
    F(10) = C0 - C0_current;
    F(11) = epsilon_cf + betaf * (Rho - Rho0) - alphaf * (Sigma_tf);
    F(12) = Sigma_tf - (Rho + E_cf * epsilon_cf);
    F(13) = Sigma - Sigma_tf;
    F(14) = Sigma_col - E_col0 * epsilon_col - ((E_col0 * lcol) / (mcol + 1)) * epsilon_col ^ (mcol + 1);
    F(15) = Sigma_tf - E_t0f * epsilon_tf - ((E_t0f * lf) / (mf + 1)) * epsilon_tf ^ (mf + 1);
end
