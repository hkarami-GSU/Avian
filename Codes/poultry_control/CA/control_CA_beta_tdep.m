clc; clear; close all;
format long
warning('off', 'all')
tic;

if ~exist('Results', 'dir')
mkdir('Results');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Needs modification (CA or MN)
m = 25; 

theta0_default = zeros(3*m, 1);
theta0_default(1) = 1-1e-1;
theta0_default(m+1) = 1-1e-1;
theta0_default(2*m+1) = 1-1e-1;

lambda_vec_1 = [0.04,0.11,0.18,0.25,0.32,0.45,0.6,0.9,1.1];
lambda_vec_2 = [0.02,0.06,0.11,0.15,0.2,0.4,0.6,0.8,1];
lambda_vec_3 = [0.24,0.72,1.25,1.8,2.4,2.7,4.5,6.5,8];

Np = [5575; 621; 114]; 

beta_data = readtable('transmission_rates.xlsx');
t_beta = double(beta_data{:, 1}); 
beta_group1 = double(beta_data{:, 2});  
beta_group2 = double(beta_data{:, 3}); 
beta_group3 = double(beta_data{:, 4});  
beta_data_groups = {t_beta, beta_group1; t_beta, beta_group2; t_beta, beta_group3};

alpha11 = 2.40e-2;
alpha12 = 6.98e-2;
alpha13 = 1.20e-3;
alpha22 = 5.25e-3;
alpha23 = 4.56e-3;
alpha33 = 5.50e-1;
alpha_params = [alpha11; alpha12; alpha13; alpha22; alpha23; alpha33];

n = 3; 
alpha_matrix = build_alpha_matrix(alpha_params, n);
Iw_data_all = readtable('wildbird_Iw_CA_all.xlsx');
I0_p = [1; 1; 7];
options1 = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt',...
'InitDamping',1e2,'Display','iter','TolX',1e-12,'TolFun',1e-12,...
'MaxIterations',20);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
a = 1;
b = 36;  
n = 3;  
phi = 30;
delta = 30/28;  
delta_vec = [delta; delta; delta];
step = 1;
tau = (a:step:b)';

tw = Iw_data_all.Time;
if isdatetime(tw) || isduration(tw)
tw = (1:length(tw))';
end
tw = double(tw(:));

num_iter = 100;
I_mat = zeros(length(Iw_data_all.Time), num_iter);
for k = 1:num_iter
field_name = sprintf('I_iter_%d', k);
I_mat(:, k) = double(Iw_data_all.(field_name)(:)); 
end
I_w = mean(I_mat, 2);  

obs_data = readtable('observation_data.xlsx');
obs_incidence = [obs_data.Poultry_Group1, obs_data.Poultry_Group2, obs_data.Poultry_Group3];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial conditions
S0_p = Np - I0_p;
C0_p = zeros(n,1);
Y0 = [S0_p; I0_p; C0_p];
Pn = [-ones(n,1); zeros(n,1); zeros(n,1)];  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Colors for plotting
colors = [
0.55, 0.00, 0.55;   % Dark Magenta
0.85, 0.33, 0.10;   % Dark Orange
0.87, 0.72, 0.53;   % Burlywood
];

legend_labels = {
'Control 1', ...
'Control 2', ...
'Control 3'
};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_all_all = cell(length(lambda_vec_1), 1);
I_all_all = cell(length(lambda_vec_1), 1);
C_all_all = cell(length(lambda_vec_1), 1);
u_all_all = cell(length(lambda_vec_1), 1);

num_lambdas = length(lambda_vec_1);
num_controls = 3;
total_runs = num_lambdas * num_controls;

lambda_all = cell(total_runs, 1);
control_idx_all = zeros(total_runs, 1);
for idx = 1:total_runs
li = ceil(idx / num_controls);
lambda_all{idx} = [lambda_vec_1(li); lambda_vec_2(li); lambda_vec_3(li)];
control_idx_all(idx) = mod(idx - 1, num_controls) + 1;
end

S_flat = cell(total_runs, 1);
I_flat = cell(total_runs, 1);
C_flat = cell(total_runs, 1);
u_flat = cell(total_runs, 1);

parfor idx = 1:total_runs
lambda = lambda_all{idx};
num_control_fun = control_idx_all(idx);
theta1 = theta0_default;

[theta, ~, ~, ~, ~, ~, ~] = lsqnonlin(@(theta) U(theta, beta_data_groups, ...
alpha_matrix, phi, delta_vec, lambda, m, tau, Y0, Pn, a, b, Np, ...
I_w, tw, n, num_control_fun), theta1, [], [], options1);

[~, Y] = ode15s(@(t,y) farm_network(t, y, theta, beta_data_groups, n, phi, ...
delta_vec, Np, alpha_matrix, I_w, tw, a, b, m), tau, Y0);

S_flat{idx} = Y(:, 1:n);
I_flat{idx} = Y(:, n+1:2*n);
C_flat{idx} = Y(:, 2*n+1:3*n);

u_vals = zeros(length(tau), n);
for k = 1:n
u_vals(:,k) = u(tau, theta((k-1)*m+1:k*m), m, a, b);
end
u_flat{idx} = u_vals;
end

for lambda_idx = 1:num_lambdas
S_temp = cell(3, 1);
I_temp = cell(3, 1);
C_temp = cell(3, 1);
u_temp = cell(3, 1);
for num_control_fun = 1:num_controls
idx = (lambda_idx - 1) * num_controls + num_control_fun;
S_temp{num_control_fun} = S_flat{idx};
I_temp{num_control_fun} = I_flat{idx};
C_temp{num_control_fun} = C_flat{idx};
u_temp{num_control_fun} = u_flat{idx};
end
S_all_all{lambda_idx} = S_temp;
I_all_all{lambda_idx} = I_temp;
C_all_all{lambda_idx} = C_temp;
u_all_all{lambda_idx} = u_temp;
end

fprintf('Parallel optimization complete for all %d lambda values\n', length(lambda_vec_1));
fprintf('Total runtime: %.2f seconds\n', toc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_date = datetime(2022, 8, 1);
tau_dates = start_date + calmonths(tau - 1);
t_start = tau_dates(1);
t_end = tau_dates(end);
t_mid = t_start + (t_end - t_start)/2;
xticklabel_fontsize = 13;
markers = {'o', 's', 'd', '^', 'p'};
Iw_interp = interp1(tw, I_w, tau, 'spline');

Incidence_all_all = cell(num_lambdas, 1);
Cumulative_all_all = cell(num_lambdas, 1);
J_all_all = cell(num_lambdas, 1);

for lambda_idx = 1:num_lambdas
lambda = [lambda_vec_1(lambda_idx); lambda_vec_2(lambda_idx); lambda_vec_3(lambda_idx)];
S_all = S_all_all{lambda_idx};
I_all = I_all_all{lambda_idx};
C_all = C_all_all{lambda_idx};
u_all = u_all_all{lambda_idx};

Incidence_all = cell(3, 1);
Cumulative_all = cell(3, 1);

for scen = 1:3
if isempty(I_all{scen}) || isempty(u_all{scen})
continue;
end

S_curr = S_all{scen};
I_curr = I_all{scen};
u_curr = u_all{scen};

inc_curr = zeros(length(tau), n);
for k_t = 1:length(tau)
for i = 1:n
beta_t = spline(beta_data_groups{i,1}, beta_data_groups{i,2}, tau(k_t));
spatial_inf = 0;
for j = 1:n
spatial_inf = spatial_inf + (I_curr(k_t, j) / Np(j)) * alpha_matrix(i, j);
end
inc_curr(k_t, i) = (beta_t * (1 - u_curr(k_t, i)) * Iw_interp(k_t) + spatial_inf) * S_curr(k_t, i);
end
end
Incidence_all{scen} = inc_curr;

cum_curr = zeros(length(tau), n);
for i = 1:n
cum_curr(:, i) = cumsum(inc_curr(:, i));
end
Cumulative_all{scen} = cum_curr;
end

J_all = zeros(3, 1);
for scen = 1:3
if isempty(u_all{scen}) || isempty(S_all{scen}) || isempty(C_all{scen})
continue;
end
u_curr = u_all{scen};
S_curr = S_all{scen};
C_curr = C_all{scen};

c_vals = zeros(length(tau), n);
for k = 1:n
u_k = u_curr(:, k);
switch scen
case 1
c_vals(:, k) = -u_k - log(1 - u_k);
case 2
c_vals(:, k) = -log(1 - u_k.^2);
case 3
c_vals(:, k) = -u_k .* log(1 - u_k);
end
end

J_scen = 0;
for k = 1:n
S_diff = S_curr(1, k) - S_curr(end, k);
control_cost = trapz(tau, lambda(k) * c_vals(:, k));
carrier_cost = trapz(tau, delta * C_curr(:, k));
J_scen = J_scen + S_diff + control_cost + carrier_cost;
end
J_all(scen) = J_scen;
end

Incidence_all_all{lambda_idx} = Incidence_all;
Cumulative_all_all{lambda_idx} = Cumulative_all;
J_all_all{lambda_idx} = J_all;

fprintf('\n');
fprintf('===== Lambda Set %d: [%.2f, %.2f, %.2f] =====\n', lambda_idx, lambda(1), lambda(2), lambda(3));
fprintf('  Control        |   Cost  \n');
fprintf('-------------------------------------\n');
for scen = 1:3
fprintf('  Control Cost %d |   %12.6f    \n', scen, J_all(scen));
end
fprintf('=====================================\n');

mat_filename = sprintf('Results/Poultry_Control_Results_tdep_%.2f_%.2f_%.2f.mat', lambda(1), lambda(2), lambda(3));
save(mat_filename, 'tau', 'lambda', 'S_all', 'I_all', 'C_all', 'u_all', ...
'Incidence_all', 'Cumulative_all', 'obs_incidence', 'colors', 'legend_labels', 'start_date', 'J_all');
fprintf('Results saved to %s\n', mat_filename);
end

Lambda_Set = cell(num_lambdas * 3, 1);
Control = cell(num_lambdas * 3, 1);
Cost = zeros(num_lambdas * 3, 1);
row = 1;
for lambda_idx = 1:num_lambdas
lambda = [lambda_vec_1(lambda_idx); lambda_vec_2(lambda_idx); lambda_vec_3(lambda_idx)];
for scen = 1:3
Lambda_Set{row} = sprintf('[%.2f, %.2f, %.2f]', lambda(1), lambda(2), lambda(3));
Control{row} = sprintf('Control %d', scen);
Cost(row) = J_all_all{lambda_idx}(scen);
row = row + 1;
end
end
Cost_Table = table(Lambda_Set, Control, Cost);
writetable(Cost_Table, 'Results/Cost_Functional_Table.xlsx');
fprintf('Cost table saved to Results/Cost_Functional_Table.xlsx\n');

lambda_idx = 1;
lambda = [lambda_vec_1(lambda_idx); lambda_vec_2(lambda_idx); lambda_vec_3(lambda_idx)];
S_all = S_all_all{lambda_idx};
I_all = I_all_all{lambda_idx};
C_all = C_all_all{lambda_idx};
u_all = u_all_all{lambda_idx};
Incidence_all = Incidence_all_all{lambda_idx};
Cumulative_all = Cumulative_all_all{lambda_idx};
J_all = J_all_all{lambda_idx};

figure(1);
set(gcf, 'Position', [100, 100, 1200, 800]);

for k = 1:n
subplot(3, 3, k)
hold on;
for scen = 1:length(S_all)
if isempty(S_all{scen}), continue; end
plot(tau_dates, S_all{scen}(:, k), 'LineWidth', 2, 'Color', colors(scen, :), ...
'DisplayName', legend_labels{scen});
end
title(sprintf('$S_%d$', k), 'Interpreter', 'latex');
set(gca, 'fontsize', 18, 'TickLabelInterpreter', 'latex');
xlim([t_start t_end]);
xticks([t_start t_mid t_end]);
xticklabels([]);
hold off;

subplot(3, 3, k+3)
hold on;
for scen = 1:length(I_all)
if isempty(I_all{scen}), continue; end
plot(tau_dates, I_all{scen}(:, k), 'LineWidth', 2, 'Color', colors(scen, :), ...
'DisplayName', legend_labels{scen});
end
title(sprintf('$I_%d$', k), 'Interpreter', 'latex');
set(gca, 'fontsize', 18, 'TickLabelInterpreter', 'latex');
if k == 1
legend('show', 'Location', 'northeast', 'Interpreter', 'latex', 'Color', 'none');
end
xlim([t_start t_end]);
xticks([t_start t_mid t_end]);
xticklabels([]);
hold off;

subplot(3, 3, k+6)
hold on;
for scen = 1:length(C_all)
if isempty(C_all{scen}), continue; end
plot(tau_dates, C_all{scen}(:, k), 'LineWidth', 2, 'Color', colors(scen, :), ...
'DisplayName', legend_labels{scen});
end
title(sprintf('$C_%d$', k), 'Interpreter', 'latex');
set(gca, 'fontsize', 18, 'TickLabelInterpreter', 'latex');
xlim([t_start t_end]);
xticks([t_start t_mid t_end]);
xticklabels([]);
yl = ylim;
text(t_start, yl(1), string(t_start, 'MMM yyyy'), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', xticklabel_fontsize);
text(t_mid, yl(1), string(t_mid, 'MMM yyyy'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', xticklabel_fontsize);
text(t_end, yl(1), string(t_end, 'MMM yyyy'), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', xticklabel_fontsize);
hold off;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2: Controls (u)
figure(2);
set(gcf, 'Position', [150, 150, 1200, 400]); 

for k = 1:n
subplot(1, n, k)
hold on;
for scen = 1:length(u_all)
if isempty(u_all{scen}), continue; end
plot(tau_dates, u_all{scen}(:, k), 'LineWidth', 2, 'Color', colors(scen, :), ...
'DisplayName', legend_labels{scen});
end
title(sprintf('$u_%d$', k), 'Interpreter', 'latex');
set(gca, 'fontsize', 18, 'TickLabelInterpreter', 'latex');
if k == 1
legend('show', 'Location', 'northeast', 'Interpreter', 'latex', 'Color', 'none');
end
xlim([t_start t_end]);
ylim([0, 1]);
xticks([t_start t_mid t_end]);
xticklabels([]);
yl = ylim;
text(t_start, yl(1), string(t_start, 'MMM yyyy'), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', xticklabel_fontsize);
text(t_mid, yl(1), string(t_mid, 'MMM yyyy'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', xticklabel_fontsize);
text(t_end, yl(1), string(t_end, 'MMM yyyy'), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', xticklabel_fontsize);
hold off;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3: Incidence (curves)
figure(3);
set(gcf, 'Position', [200, 200, 1200, 400]);

for k = 1:n
subplot(1, n, k);
hold on;

len_obs = min(length(tau), size(obs_incidence, 1));
plot(tau_dates(1:len_obs), obs_incidence(1:len_obs, k), 'ko', 'MarkerSize', 6, 'LineWidth', 2, 'DisplayName', 'Reported Data');

for scen = 1:length(Incidence_all)
if isempty(Incidence_all{scen}), continue; end
plot(tau_dates, Incidence_all{scen}(:, k), 'LineWidth', 2, 'Color', colors(scen, :), ...
'DisplayName', legend_labels{scen});
end

set(gca, 'fontsize', 18, 'TickLabelInterpreter', 'latex');
if k == 1
ylabel('Incidence', 'Interpreter', 'latex', 'FontSize', 18);
legend('show', 'Location', 'northeast', 'Interpreter', 'latex', 'Color', 'none');
end
xlim([t_start t_end]);
xticks([t_start t_mid t_end]);
xticklabels([]);
yl = ylim;
text(t_start, yl(1), string(t_start, 'MMM yyyy'), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', xticklabel_fontsize);
text(t_mid, yl(1), string(t_mid, 'MMM yyyy'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', xticklabel_fontsize);
text(t_end, yl(1), string(t_end, 'MMM yyyy'), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', xticklabel_fontsize);
hold off;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4: Incidence Integer (circles at observation points)
figure(4);
set(gcf, 'Position', [250, 250, 1200, 400]);

for k = 1:n
subplot(1, n, k);
hold on;

len_obs = min(length(tau), size(obs_incidence, 1));

num_scen = length(Incidence_all);
offsets = linspace(-5, 5, num_scen + 1);

plot(tau_dates(1:len_obs), obs_incidence(1:len_obs, k), 'ko', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Reported Data');

for scen = 1:length(Incidence_all)
if isempty(Incidence_all{scen}), continue; end
inc_rounded = round(Incidence_all{scen}(:, k));
marker_style = markers{mod(scen-1, length(markers)) + 1};
x_offset = tau_dates(1:len_obs) + days(offsets(scen));
plot(x_offset, inc_rounded(1:len_obs), marker_style, 'Color', colors(scen, :), ...
'MarkerFaceColor', colors(scen, :), 'MarkerSize', 6, 'DisplayName', legend_labels{scen});
end

set(gca, 'fontsize', 18, 'TickLabelInterpreter', 'latex');
if k == 1
ylabel('Incidence', 'Interpreter', 'latex', 'FontSize', 18);
legend('show', 'Location', 'northeast', 'Interpreter', 'latex', 'Color', 'none');
end
xlim([t_start t_end]);
xticks([t_start t_mid t_end]);
xticklabels([]);
yl = ylim;
text(t_start, yl(1), string(t_start, 'MMM yyyy'), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', xticklabel_fontsize);
text(t_mid, yl(1), string(t_mid, 'MMM yyyy'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', xticklabel_fontsize);
text(t_end, yl(1), string(t_end, 'MMM yyyy'), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', xticklabel_fontsize);
hold off;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5: Cumulative Incidence (curves)
figure(5);
set(gcf, 'Position', [300, 300, 1200, 400]);

for k = 1:n
subplot(1, n, k);
hold on;

cum_obs = cumsum(obs_incidence(:, k));
len_obs = min(length(tau), size(obs_incidence, 1));
plot(tau_dates(1:len_obs), cum_obs(1:len_obs), 'ko', 'MarkerSize', 6, 'DisplayName', 'Reported Data');

for scen = 1:length(Cumulative_all)
if isempty(Cumulative_all{scen}), continue; end
plot(tau_dates, Cumulative_all{scen}(:, k), 'LineWidth', 2, 'Color', colors(scen, :), ...
'DisplayName', legend_labels{scen});
end

set(gca, 'fontsize', 18, 'TickLabelInterpreter', 'latex');
if k == 1
ylabel('Cumulative Cases', 'Interpreter', 'latex', 'FontSize', 18);
legend('show', 'Location', 'best', 'Interpreter', 'latex', 'Color', 'none');
end
xlim([t_start t_end]);
xticks([t_start t_mid t_end]);
xticklabels([]);
yl = ylim;
text(t_start, yl(1), string(t_start, 'MMM yyyy'), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', xticklabel_fontsize);
text(t_mid, yl(1), string(t_mid, 'MMM yyyy'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', xticklabel_fontsize);
text(t_end, yl(1), string(t_end, 'MMM yyyy'), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', xticklabel_fontsize);
hold off;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6: Cumulative Incidence Integer (circles)
figure(6);
set(gcf, 'Position', [350, 350, 1200, 400]);

for k = 1:n
subplot(1, n, k);
hold on;

cum_obs = cumsum(obs_incidence(:, k));
len_obs = min(length(tau), size(obs_incidence, 1));

% Horizontal offset in days for each scenario to avoid overlap
num_scen = length(Cumulative_all);
offsets = linspace(-5, 5, num_scen + 1);

plot(tau_dates(1:len_obs), cum_obs(1:len_obs), 'ko', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Reported Data');

for scen = 1:length(Cumulative_all)
if isempty(Cumulative_all{scen}), continue; end
cum_rounded = round(Cumulative_all{scen}(:, k));
marker_style = markers{mod(scen-1, length(markers)) + 1};
x_offset = tau_dates(1:len_obs) + days(offsets(scen));
plot(x_offset, cum_rounded(1:len_obs), marker_style, 'Color', colors(scen, :), ...
'MarkerFaceColor', colors(scen, :), 'MarkerSize', 6, 'DisplayName', legend_labels{scen});
end

set(gca, 'fontsize', 18, 'TickLabelInterpreter', 'latex');
if k == 1
ylabel('Cumulative Cases', 'Interpreter', 'latex', 'FontSize', 18);
end
if k == 1
lgd = legend('show', 'Location', 'southeast', 'Interpreter', 'latex', 'Color', 'none');
pos = lgd.Position;
pos(1) = pos(1) + 0.01;
lgd.Position = pos;
end
xlim([t_start t_end]);
xticks([t_start t_mid t_end]);
xticklabels([]);
yl = ylim;
text(t_start, yl(1), string(t_start, 'MMM yyyy'), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', xticklabel_fontsize);
text(t_mid, yl(1), string(t_mid, 'MMM yyyy'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', xticklabel_fontsize);
text(t_end, yl(1), string(t_end, 'MMM yyyy'), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', xticklabel_fontsize);
hold off;
end

fprintf('Figures generated.\n');
fprintf('All %d lambda combinations processed and saved.\n', num_lambdas);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION DEFINITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = farm_network(t, y, theta, beta_data_groups, n, phi, delta, Np, alpha_matrix, I_w, tw, a, b, m)
S = y(1:n);
I = y(n+1:2*n);
C = y(2*n+1:3*n);
I_w_t = interp1(tw, I_w, t, 'spline');

u_vals = zeros(n,1);
for k = 1:n
u_vals(k)=u(t,theta((k-1)*m+1:k*m),m,a,b);
end

dSdt = zeros(n,1);
dIdt = zeros(n,1);
dCdt = zeros(n,1);

for i = 1:n
beta_t = spline(beta_data_groups{i,1}, beta_data_groups{i,2}, t);
spatial_infection = 0;
for j = 1:n
spatial_infection = spatial_infection + (I(j)/Np(j)) * alpha_matrix(i,j);
end
dSdt(i)=-beta_t*(1-u_vals(i))*S(i)*I_w_t-S(i)*spatial_infection+delta(i)*C(i);
dIdt(i)=beta_t*(1-u_vals(i))*S(i)*I_w_t+S(i)*spatial_infection-phi*I(i);
dCdt(i)=phi*I(i)-delta(i)*C(i);
end
dydt = [dSdt; dIdt; dCdt];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dpdt = p_adj(t, p, beta_data_groups, alpha_matrix, phi, delta_vec, theta, m, tau, a, b, S, I, Np, I_w, tw, n)
SS = zeros(n,1);
II = zeros(n,1);
for i = 1:n
SS(i) = spline(tau, S(:,i), t);
II(i) = spline(tau, I(:,i), t);
end
Iw_t = interp1(tw, I_w, t, 'spline');
dpdt = zeros(3*n, 1);

u_vals = zeros(n,1);
for k = 1:n
u_vals(k) = u(t, theta((k-1)*m+1:k*m), m, a, b);
end

for k = 1:n
beta_k = spline(beta_data_groups{k,1}, beta_data_groups{k,2}, t);

p1_k = p((k-1)*3+1);
p2_k = p((k-1)*3+2);
p3_k = p((k-1)*3+3);

% dp1^(k)/dt
sum_term = 0;
for j = 1:n
sum_term = sum_term+alpha_matrix(k,j)*(II(j)/Np(j));
end
dpdt((k-1)*3+1)=(beta_k*(1-u_vals(k))*Iw_t+sum_term)*(p1_k-p2_k);

% dp2^(k)/dt
sum_term2 = 0;
for i = 1:n
p1_i = p((i-1)*3 + 1);
p2_i = p((i-1)*3 + 2);
sum_term2 = sum_term2 + alpha_matrix(i,k) * (SS(i)/Np(k)) * (p1_i - p2_i);
end
dpdt((k-1)*3 + 2) = sum_term2 + phi * (p2_k - p3_k);

% dp3^(k)/dt
dpdt((k-1)*3 + 3) = delta_vec(k)*(p3_k - p1_k - 1);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q = U(theta, beta_data_groups, alpha_matrix, phi, delta_vec, lambda, m, tau, Y0, Pn, a, b, Np, I_w, tw, n, num_control_fun)
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[~,Y] = ode15s(@(t,y) farm_network(t,y,theta,beta_data_groups,n,phi,delta_vec,Np,alpha_matrix,...
I_w, tw, a, b, m), tau, Y0, options);
S = Y(:,1:n);
I = Y(:,n+1:2*n);

[~,P] = ode15s(@(t,p) p_adj(t, p, beta_data_groups, alpha_matrix, phi, delta_vec, theta,...
m, tau, a, b, S, I, Np, I_w, tw, n), flip(tau), Pn, options);
P = flip(P);

q = zeros(length(tau), n);
Iw_vals = interp1(tw, I_w, tau, 'spline');

for k = 1:n
p1_k = P(:,(k-1)*3+1);
p2_k = P(:,(k-1)*3+2);

beta_k_vals = spline(beta_data_groups{k,1}, beta_data_groups{k,2}, tau);
beta_k_vals = beta_k_vals(:);

q(:, k) = lambda(k) .* g(tau, theta((k-1)*m+1:k*m), m, a, b, num_control_fun) + ...
beta_k_vals .* S(:,k) .* Iw_vals .* (p1_k - p2_k);
end

q = q(:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function yy = u(t, theta, m, a, b)
yy = 0;
for j = 1:m
yy = yy + theta(j) * leg(j-1, t, a, b);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = leg(j, t, a, b)
x = (2.*t - a - b)./(b - a);

if j == 0
P = ones(size(t));
elseif j == 1
P = x;
else
P1 = ones(size(t)); P2 = x;
for k = 2:j
P3 = ((2*(k-1)+1).*x.*P2 - (k-1).*P1)./k;
P1 = P2; P2 = P3;
end
P = P3;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Control = g(tau, theta, m, a, b, num_control_fun)
switch num_control_fun
case 1  % c'(u) = u/(1-u)
u_vals = u(tau, theta, m, a, b);
Control = u_vals ./ (1 - u_vals);
case 2  % c'(u) = 2u/(1-u^2)
u_vals = u(tau, theta, m, a, b);
Control = 2 * u_vals ./ (1 - u_vals.^2);
case 3  % c'(u) = -ln(1-u) + u/(1-u)
u_vals = u(tau, theta, m, a, b);
Control = -log(1 - u_vals) + u_vals ./ (1 - u_vals);
otherwise
error('Invalid control function number!');
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alpha_matrix = build_alpha_matrix(alpha, n)
alpha_matrix = zeros(n, n);
idx = 1;
for i = 1:n
for j = i:n
alpha_matrix(i,j) = alpha(idx);
alpha_matrix(j,i) = alpha(idx);
idx = idx + 1;
end
end
end
