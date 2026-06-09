clc; clear; close all;
format long
warning('off', 'all')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration parameters
m_beta = 10; % Number of Legendre polynomials for each beta_i(t)
noise_level_incidence = 0.05;    % 5% noise
noise_level_cumulative = 0.05;   % 5% noise
t0 = 5*pi/6;
n = 3;
delta_min = 30/60;
delta_max = 30/21;
phi = 30;
Np = [5575; 621; 114];
NumCurves = 100;
MaxCurves = 100;
cumulative_downweight = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load observation data
Data = readtable('observation_data.xlsx');
tdata = Data.Months;
Data_p1_incidence = Data.Poultry_Group1;
Data_p2_incidence = Data.Poultry_Group2;
Data_p3_incidence = Data.Poultry_Group3;

if isdatetime(tdata) || isduration(tdata)
tdata = (1:length(tdata))';
end
tdata = double(tdata(:));
Data_p1_cumulative = cumsum(Data_p1_incidence);
Data_p2_cumulative = cumsum(Data_p2_incidence);
Data_p3_cumulative = cumsum(Data_p3_incidence);

a = min(tdata);
b = max(tdata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalization 
mean_p1_incidence = mean(abs(Data_p1_incidence(Data_p1_incidence > 0)));
mean_p2_incidence = mean(abs(Data_p2_incidence(Data_p2_incidence > 0)));
mean_p3_incidence = mean(abs(Data_p3_incidence(Data_p3_incidence > 0)));
mean_p1_cumulative = mean(abs(Data_p1_cumulative(Data_p1_cumulative > 0)));
mean_p2_cumulative = mean(abs(Data_p2_cumulative(Data_p2_cumulative > 0)));
mean_p3_cumulative = mean(abs(Data_p3_cumulative(Data_p3_cumulative > 0)));

% Calculate absolute noise standard deviations based on relative noise levels
noise_p1_incidence = noise_level_incidence * mean_p1_incidence;
noise_p2_incidence = noise_level_incidence * mean_p2_incidence;
noise_p3_incidence = noise_level_incidence * mean_p3_incidence;
noise_p1_cumulative = noise_level_cumulative * mean_p1_cumulative;
noise_p2_cumulative = noise_level_cumulative * mean_p2_cumulative;
noise_p3_cumulative = noise_level_cumulative * mean_p3_cumulative;

% Calculate normalization weights (inverse of mean)
weight_incidence_1 = 1 / mean_p1_incidence;
weight_cumulative_1 = 1 / mean_p1_cumulative;
weight_incidence_2 = 1 / mean_p2_incidence;
weight_cumulative_2 = 1 / mean_p2_cumulative;
weight_incidence_3 = 1 / mean_p3_incidence;
weight_cumulative_3 = 1 / mean_p3_cumulative;

% Downweight cumulative data 
weight_cumulative_1 = weight_cumulative_1 * cumulative_downweight;
weight_cumulative_2 = weight_cumulative_2 * cumulative_downweight;
weight_cumulative_3 = weight_cumulative_3 * cumulative_downweight;

% Load wildbird data with all iterations
Iw_data_all = readtable('wildbird_Iw_CA_all.xlsx');
tw = Iw_data_all.Time;

if isdatetime(tw) || isduration(tw)
    tw = (1:length(tw))';
end
tw = double(tw(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial conditions
I0_p = [Data_p1_incidence(1); Data_p2_incidence(1); Data_p3_incidence(1)];
S0_p = Np - I0_p;
X0 = [S0_p; I0_p; zeros(n,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate noisy observation curves
ndata = length(tdata);

curves_p1_incidence = zeros(ndata, MaxCurves);
curves_p2_incidence = zeros(ndata, MaxCurves);
curves_p3_incidence = zeros(ndata, MaxCurves);
curves_p1_cumulative = zeros(ndata, MaxCurves);
curves_p2_cumulative = zeros(ndata, MaxCurves);
curves_p3_cumulative = zeros(ndata, MaxCurves);

rng(0);

for i = 1:MaxCurves
noisy_p1_incidence = Data_p1_incidence + normrnd(0, noise_p1_incidence, ndata, 1);
noisy_p1_incidence(noisy_p1_incidence < 0) = 0;
noisy_p1_cumulative = Data_p1_cumulative + normrnd(0, noise_p1_cumulative, ndata, 1);
noisy_p1_cumulative(noisy_p1_cumulative < 0) = 0;
curves_p1_incidence(:, i) = noisy_p1_incidence;
curves_p1_cumulative(:, i) = noisy_p1_cumulative;

noisy_p2_incidence = Data_p2_incidence + normrnd(0, noise_p2_incidence, ndata, 1);
noisy_p2_incidence(noisy_p2_incidence < 0) = 0;
noisy_p2_cumulative = Data_p2_cumulative + normrnd(0, noise_p2_cumulative, ndata, 1);
noisy_p2_cumulative(noisy_p2_cumulative < 0) = 0;
curves_p2_incidence(:, i) = noisy_p2_incidence;
curves_p2_cumulative(:, i) = noisy_p2_cumulative;

noisy_p3_incidence = Data_p3_incidence + normrnd(0, noise_p3_incidence, ndata, 1);
noisy_p3_incidence(noisy_p3_incidence < 0) = 0;
noisy_p3_cumulative = Data_p3_cumulative + normrnd(0, noise_p3_cumulative, ndata, 1);
noisy_p3_cumulative(noisy_p3_cumulative < 0) = 0;
curves_p3_incidence(:, i) = noisy_p3_incidence;
curves_p3_cumulative(:, i) = noisy_p3_cumulative;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial guesses and bounds
beta1_0 = [0.0001, 0, 0, 0, 0, 0, 0, 0, 0, 0];
beta2_0 = [0.0001, 0, 0, 0, 0, 0, 0, 0, 0, 0];
beta3_0 = [0.0001, 0, 0, 0, 0, 0, 0, 0, 0, 0];
alpha0 = [0.05, 0.01, 0.02, 0.02, 0.01, 0.5];

theta0 = [beta1_0; beta2_0; beta3_0; alpha0];

lb = [-Inf*ones(n*m_beta,1); zeros(6,1)];        
ub = [Inf*ones(n*m_beta,1); Inf*ones(6,1)];      

options = optimoptions('lsqcurvefit', ...
    'Algorithm', 'levenberg-marquardt', ...
    'Display', 'iter', ...
    'MaxIterations', 1000, ...
    'TolX', 1e-7, ...
    'TolFun', 1e-7, ...
    'InitDamping',1e10,...
    'UseParallel', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN
Beta_params_temp = zeros(MaxCurves, n*m_beta);
Alpha_params_temp = zeros(MaxCurves, 6);
Delta_params_temp = zeros(MaxCurves, 1);
curves_results_p1_inc_temp = zeros(MaxCurves, ndata);
curves_results_p2_inc_temp = zeros(MaxCurves, ndata);
curves_results_p3_inc_temp = zeros(MaxCurves, ndata);
curves_results_p1_cum_temp = zeros(MaxCurves, ndata);
curves_results_p2_cum_temp = zeros(MaxCurves, ndata);
curves_results_p3_cum_temp = zeros(MaxCurves, ndata);
success_flags_p = false(MaxCurves, 1);
ExpData = zeros(2*n*length(tdata), MaxCurves);
Delta_values = zeros(MaxCurves, 1);
Iw_columns = cell(MaxCurves, 1);

parfor i = 1:MaxCurves
ExpData(:,i) = [curves_p1_incidence(:,i)*weight_incidence_1; curves_p1_cumulative(:,i)*weight_cumulative_1; ...
curves_p2_incidence(:,i)*weight_incidence_2; curves_p2_cumulative(:,i)*weight_cumulative_2; ...
curves_p3_incidence(:,i)*weight_incidence_3; curves_p3_cumulative(:,i)*weight_cumulative_3];
end

rng(0);
for i = 1:MaxCurves
Delta_values(i) = delta_min + (delta_max - delta_min) * rand();
column_name = sprintf('I_iter_%d', i);
Iw_columns{i} = double(Iw_data_all.(column_name)(:));
end

parfor i = 1:MaxCurves
ObsData_p = ExpData(:,i);
delta_i = Delta_values(i);
I_w_i = Iw_columns{i};

theta_opt = lsqcurvefit(@(theta, tdata) Phi_poultry(theta, tdata, X0, n, m_beta, ...
phi, delta_i, Np, I_w_i, tw, a, b, weight_incidence_1, weight_cumulative_1, ...
weight_incidence_2, weight_cumulative_2, weight_incidence_3, weight_cumulative_3), ...
theta0, tdata, ObsData_p, lb, ub, options);

Beta_params_temp(i,:) = theta_opt(1:n*m_beta)';
Alpha_params_temp(i,:) = theta_opt(n*m_beta+1:n*m_beta+6)';
Delta_params_temp(i) = delta_i;

q_all = I_poultry(theta_opt, tdata, X0, n, m_beta, phi, delta_i, Np, I_w_i, tw, a, b);
curves_results_p1_inc_temp(i,:) = q_all(1:ndata)';
curves_results_p1_cum_temp(i,:) = q_all(ndata+1:2*ndata)';
curves_results_p2_inc_temp(i,:) = q_all(2*ndata+1:3*ndata)';
curves_results_p2_cum_temp(i,:) = q_all(3*ndata+1:4*ndata)';
curves_results_p3_inc_temp(i,:) = q_all(4*ndata+1:5*ndata)';
curves_results_p3_cum_temp(i,:) = q_all(5*ndata+1:end)';
success_flags_p(i) = true;
end

success_indices_p = find(success_flags_p);
if length(success_indices_p) < NumCurves
error('Insufficient successful fits: %d/%d', length(success_indices_p), NumCurves);
end

Beta_params = Beta_params_temp(success_indices_p(1:NumCurves), :);
Alpha_params = Alpha_params_temp(success_indices_p(1:NumCurves), :);
Delta_params = Delta_params_temp(success_indices_p(1:NumCurves));
curves_results_p1_inc = curves_results_p1_inc_temp(success_indices_p(1:NumCurves), :);
curves_results_p1_cum = curves_results_p1_cum_temp(success_indices_p(1:NumCurves), :);
curves_results_p2_inc = curves_results_p2_inc_temp(success_indices_p(1:NumCurves), :);
curves_results_p2_cum = curves_results_p2_cum_temp(success_indices_p(1:NumCurves), :);
curves_results_p3_inc = curves_results_p3_inc_temp(success_indices_p(1:NumCurves), :);
curves_results_p3_cum = curves_results_p3_cum_temp(success_indices_p(1:NumCurves), :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta_params_opt = mean(Beta_params, 1)';
beta_params_opt = reshape(beta_params_opt, m_beta, n)'; % Each row is one group's coefficients
alpha_params_opt = mean(Alpha_params, 1)';
delta_opt = mean(Delta_params);

alpha_opt = zeros(n, n);
idx = 1;
for i = 1:n
for j = i:n
alpha_opt(i,j) = alpha_params_opt(idx);
alpha_opt(j,i) = alpha_params_opt(idx);
idx = idx + 1;
end
end

theta_opt = [beta_params_opt(:); alpha_params_opt];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1: Poultry cumulative fit
figure1 = figure('Units','inches','Position',[1 1 15 5]);

start_date = datetime(2022, 8, 1);
tdata_plot = start_date + calmonths(tdata - 1);
t_start = tdata_plot(1);
t_end   = tdata_plot(end);
t_mid   = t_start + (t_end - t_start)/2;
font_date = 14;

subplot(1,3,1);
h1 = plot(tdata_plot, curves_p1_cumulative, 'Color', [0, 1, 1], 'LineWidth', 1);
hold on;
h2 = plot(tdata_plot, curves_results_p1_cum', 'Color', [0, 0.5, 0], 'LineWidth', 1);
h3 = plot(tdata_plot, mean(curves_results_p1_cum), '-', 'Color', [0 0.5 0], 'LineWidth', 2);
h4 = plot(tdata_plot, Data_p1_cumulative, 'ko', 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
ylabel('Cumulative Cases', 'Interpreter','latex','FontSize',18);
xlim([t_start t_end]);
xticks([t_start t_mid t_end]);
set(gca, 'XTickLabel', []);
text(t_start, min(ylim), datestr(t_start, 'mmm yyyy'), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
text(t_mid, min(ylim), datestr(t_mid, 'mmm yyyy'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
text(t_end, min(ylim), datestr(t_end, 'mmm yyyy'), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
legend([h1(1), h2(1), h3, h4], {'Noisy Data','Reconstructed','Mean Reconstructed','Reported Data'}, ...
    'Interpreter','latex','FontSize',12,'Location','best');
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
grid off;

subplot(1,3,2);
h1=plot(tdata_plot, curves_p2_cumulative, 'Color', [0, 1, 1], 'LineWidth', 1);
hold on;
h2=plot(tdata_plot, curves_results_p2_cum', 'Color', [1, 0.5, 0], 'LineWidth', 1);
h3=plot(tdata_plot, mean(curves_results_p2_cum), '-', 'Color', [1 0.5 0], 'LineWidth', 2);
h4=plot(tdata_plot, Data_p2_cumulative, 'ko', 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
xlim([t_start t_end]);
xticks([t_start t_mid t_end]);
set(gca, 'XTickLabel', []);
text(t_start, min(ylim), datestr(t_start, 'mmm yyyy'), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
text(t_mid, min(ylim), datestr(t_mid, 'mmm yyyy'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
text(t_end, min(ylim), datestr(t_end, 'mmm yyyy'), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
legend([h1(1), h2(1), h3, h4], {'Noisy Data','Reconstructed','Mean Reconstructed','Reported Data'}, ...
    'Interpreter','latex','FontSize',12,'Location','best');
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
grid off;

subplot(1,3,3);
h1=plot(tdata_plot, curves_p3_cumulative, 'Color', [0, 1, 1], 'LineWidth', 1);
hold on;
h2=plot(tdata_plot, curves_results_p3_cum', 'Color', [0, 0, 0.5], 'LineWidth', 1);
h3=plot(tdata_plot, mean(curves_results_p3_cum), '-', 'Color', [0 0 0.5], 'LineWidth', 2);
h4=plot(tdata_plot, Data_p3_cumulative, 'ko', 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
xlim([t_start t_end]);
xticks([t_start t_mid t_end]);
set(gca, 'XTickLabel', []);
text(t_start, min(ylim), datestr(t_start, 'mmm yyyy'), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
text(t_mid, min(ylim), datestr(t_mid, 'mmm yyyy'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
text(t_end, min(ylim), datestr(t_end, 'mmm yyyy'), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
legend([h1(1), h2(1), h3, h4], {'Noisy Data','Reconstructed','Mean Reconstructed','Reported Data'}, ...
    'Interpreter','latex','FontSize',12,'Location','best');
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
grid off;

set(figure1,'PaperUnits','inches','PaperPosition',[0 0 15 4],'PaperSize',[15 4]);
exportgraphics(figure1,'Results/poultry_cumulative_fit_independent.pdf','ContentType','vector');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2: Poultry incidence fit
figure2 = figure('Units','inches','Position',[1 1 15 5]);

subplot(1,3,1);
h1 = plot(tdata_plot, curves_p1_incidence, 'Color', [0, 1, 1], 'LineWidth', 1);
hold on;
h2 = plot(tdata_plot, curves_results_p1_inc', 'Color', [0, 0.5, 0], 'LineWidth', 1);
h3 = plot(tdata_plot, mean(curves_results_p1_inc), '-', 'Color', [0 0.5 0], 'LineWidth', 2);
h4 = plot(tdata_plot, Data_p1_incidence, 'ko', 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
ylabel('Incidence Cases', 'Interpreter','latex','FontSize',16);
xlim([t_start t_end]);
xticks([t_start t_mid t_end]);
set(gca, 'XTickLabel', []);
text(t_start, min(ylim), datestr(t_start, 'mmm yyyy'), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
text(t_mid, min(ylim), datestr(t_mid, 'mmm yyyy'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
text(t_end, min(ylim), datestr(t_end, 'mmm yyyy'), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
legend([h1(1), h2(1), h3, h4], {'Noisy Data','Reconstructed','Mean Reconstructed','Reported Data'}, ...
    'Interpreter','latex','FontSize',12,'Location','best');
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
grid off;

subplot(1,3,2);
h1=plot(tdata_plot, curves_p2_incidence, 'Color', [0, 1, 1], 'LineWidth', 1);
hold on;
h2=plot(tdata_plot, curves_results_p2_inc', 'Color', [1, 0.5, 0], 'LineWidth', 1);
h3=plot(tdata_plot, mean(curves_results_p2_inc), '-', 'Color', [1 0.5 0], 'LineWidth', 2);
h4=plot(tdata_plot, Data_p2_incidence, 'ko', 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
xlim([t_start t_end]);
xticks([t_start t_mid t_end]);
set(gca, 'XTickLabel', []);
text(t_start, min(ylim), datestr(t_start, 'mmm yyyy'), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
text(t_mid, min(ylim), datestr(t_mid, 'mmm yyyy'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
text(t_end, min(ylim), datestr(t_end, 'mmm yyyy'), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
legend([h1(1), h2(1), h3, h4], {'Noisy Data','Reconstructed','Mean Reconstructed','Reported Data'}, ...
    'Interpreter','latex','FontSize',12,'Location','best');
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
grid off;

subplot(1,3,3);
h1=plot(tdata_plot, curves_p3_incidence, 'Color', [0, 1, 1], 'LineWidth', 1);
hold on;
h2=plot(tdata_plot, curves_results_p3_inc', 'Color', [0, 0, 0.5], 'LineWidth', 1);
h3=plot(tdata_plot, mean(curves_results_p3_inc), '-', 'Color', [0 0 0.5], 'LineWidth', 2);
h4=plot(tdata_plot, Data_p3_incidence, 'ko', 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
xlim([t_start t_end]);
xticks([t_start t_mid t_end]);
set(gca, 'XTickLabel', []);
text(t_start, min(ylim), datestr(t_start, 'mmm yyyy'), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
text(t_mid, min(ylim), datestr(t_mid, 'mmm yyyy'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
text(t_end, min(ylim), datestr(t_end, 'mmm yyyy'), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
legend([h1(1), h2(1), h3, h4], {'Noisy Data','Reconstructed','Mean Reconstructed','Reported Data'}, ...
    'Interpreter','latex','FontSize',12,'Location','best');
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
grid off;

set(figure2,'PaperUnits','inches','PaperPosition',[0 0 15 4],'PaperSize',[15 4]);
exportgraphics(figure2,'Results/poultry_incidence_fit_independent.pdf','ContentType','vector');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3: Beta parameter histograms (all groups, all coefficients)
figure3 = figure('Units','inches','Position',[1 1 15 10]);
for i = 1:n
for j = 1:m_beta
idx = (i-1)*m_beta + j;
subplot(n, m_beta, idx);
h = histogram(Beta_params(:, idx), 'Normalization', 'probability');
ylims = [0, max(h.Values)*1.1];
ylim(ylims);
hold on;

mean_val = mean(Beta_params(:, idx));
ci_95 = prctile(Beta_params(:, idx), [2.5, 97.5]);

plot([mean_val, mean_val], ylim, 'r-', 'LineWidth', 2);
plot([ci_95(1), ci_95(1)], ylim, 'k--', 'LineWidth', 1.5);
plot([ci_95(2), ci_95(2)], ylim, 'k--', 'LineWidth', 1.5);

title(sprintf('$\\beta_{%d,%d} = %.2e$ \n $(%.2e, %.2e)$', i, j-1, mean_val, ci_95(1), ci_95(2)), ...
'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'FontSize', 10, 'TickLabelInterpreter', 'latex');
end
end

set(figure3,'PaperUnits','inches','PaperPosition',[0 0 15 10],'PaperSize',[15 10]);
exportgraphics(figure3,'Results/params_beta_independent.pdf','ContentType','vector');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4: Alpha parameter histograms
figure4 = figure('Units','inches','Position',[1 1 15 4]);
alpha_labels = {'\alpha_{11}', '\alpha_{12}', '\alpha_{13}', '\alpha_{22}', '\alpha_{23}', '\alpha_{33}'};
for j = 1:6
subplot(1, 6, j);
h = histogram(Alpha_params(:, j), 'Normalization', 'probability');
ylims = [0, max(h.Values)*1.1];
ylim(ylims);
hold on;

mean_val = mean(Alpha_params(:, j));
ci_95 = prctile(Alpha_params(:, j), [2.5, 97.5]);
initial_val = alpha0(j);

plot([mean_val, mean_val], ylim, 'r-', 'LineWidth', 2);
plot([ci_95(1), ci_95(1)], ylim, 'k--', 'LineWidth', 1.5);
plot([ci_95(2), ci_95(2)], ylim, 'k--', 'LineWidth', 1.5);
plot([initial_val, initial_val], ylim, 'y-', 'LineWidth', 2);

title(sprintf('$%s = %.2e$ \n $(%.2e, %.2e)$', alpha_labels{j}, mean_val, ci_95(1), ci_95(2)), ...
'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');
end

set(figure4,'PaperUnits','inches','PaperPosition',[0 0 15 4],'PaperSize',[15 4]);
exportgraphics(figure4,'Results/params_alpha.pdf','ContentType','vector');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5: Delta parameter histogram
figure5 = figure('Units','inches','Position',[1 1 5 4]);

h = histogram(Delta_params, 'Normalization', 'probability');
ylims = [0, max(h.Values)*1.1];
ylim(ylims);
hold on;

mean_val = mean(Delta_params);
ci_95 = prctile(Delta_params, [2.5, 97.5]);

plot([mean_val, mean_val], ylim, 'r-', 'LineWidth', 2);
plot([ci_95(1), ci_95(1)], ylim, 'k--', 'LineWidth', 1.5);
plot([ci_95(2), ci_95(2)], ylim, 'k--', 'LineWidth', 1.5);

title(sprintf('$\\delta = %.3f$\n $(%.3f, %.3f)$', mean_val, ci_95(1), ci_95(2)), ...
    'Interpreter', 'latex', 'FontSize', 14);
xlabel('Delta', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Probability', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');
grid on;

set(figure5,'PaperUnits','inches','PaperPosition',[0 0 5 4],'PaperSize',[5 4]);
exportgraphics(figure5,'Results/param_delta.pdf','ContentType','vector');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6: Time-dependent betas (all groups)
figure6 = figure('Units','inches','Position',[1 1 12 4]);
t_fine = linspace(a, b, 100);
t_fine_plot = start_date + days((t_fine - 1) * 30); % Convert to datetime using days

beta1_vals = zeros(length(t_fine), 1);
beta2_vals = zeros(length(t_fine), 1);
beta3_vals = zeros(length(t_fine), 1);
for k = 1:length(t_fine)
beta1_vals(k) = beta(t_fine(k), beta_params_opt(1,:), a, b);
beta2_vals(k) = beta(t_fine(k), beta_params_opt(2,:), a, b);
beta3_vals(k) = beta(t_fine(k), beta_params_opt(3,:), a, b);
end

subplot(1,3,1);
plot(t_fine_plot, beta1_vals, '-', 'Color', [0 0.5 0], 'LineWidth', 2);
ylabel('$\beta_1(t)$', 'Interpreter','latex','FontSize',18);
xlim([t_start t_end]);
xticks([t_start t_mid t_end]);
set(gca, 'XTickLabel', []);
text(t_start, min(ylim), datestr(t_start, 'mmm yyyy'), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
text(t_mid, min(ylim), datestr(t_mid, 'mmm yyyy'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
text(t_end, min(ylim), datestr(t_end, 'mmm yyyy'), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
grid on;

subplot(1,3,2);
plot(t_fine_plot, beta2_vals, '-', 'Color', [1 0.5 0], 'LineWidth', 2);
ylabel('$\beta_2(t)$', 'Interpreter','latex','FontSize',18);
xlim([t_start t_end]);
xticks([t_start t_mid t_end]);
set(gca, 'XTickLabel', []);
text(t_start, min(ylim), datestr(t_start, 'mmm yyyy'), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
text(t_mid, min(ylim), datestr(t_mid, 'mmm yyyy'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
text(t_end, min(ylim), datestr(t_end, 'mmm yyyy'), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
grid on;

subplot(1,3,3);
plot(t_fine_plot, beta3_vals, '-', 'Color', [0 0 0.5], 'LineWidth', 2);
ylabel('$\beta_3(t)$', 'Interpreter','latex','FontSize',18);
xlim([t_start t_end]);
xticks([t_start t_mid t_end]);
set(gca, 'XTickLabel', []);
text(t_start, min(ylim), datestr(t_start, 'mmm yyyy'), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
text(t_mid, min(ylim), datestr(t_mid, 'mmm yyyy'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
text(t_end, min(ylim), datestr(t_end, 'mmm yyyy'), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
grid on;

set(figure6,'PaperUnits','inches','PaperPosition',[0 0 12 4],'PaperSize',[12 4]);
exportgraphics(figure6,'Results/beta_time_dependent_independent.pdf','ContentType','vector');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save results to .mat file
save('Results/Poultry_Inverse_CA_tdep_INDEPENDENT_fix_results.mat', ...
    'tdata', 'tdata_plot', 'start_date', 't_start', 't_end', 't_mid', 'a', 'b', ...
    'Data_p1_incidence', 'Data_p2_incidence', 'Data_p3_incidence', ...
    'Data_p1_cumulative', 'Data_p2_cumulative', 'Data_p3_cumulative', ...
    'curves_p1_incidence', 'curves_p2_incidence', 'curves_p3_incidence', ...
    'curves_p1_cumulative', 'curves_p2_cumulative', 'curves_p3_cumulative', ...
    'curves_results_p1_inc', 'curves_results_p2_inc', 'curves_results_p3_inc', ...
    'curves_results_p1_cum', 'curves_results_p2_cum', 'curves_results_p3_cum', ...
    'Beta_params', 'Alpha_params', 'Delta_params', ...
    'beta_params_opt', 'alpha_params_opt', 'alpha_opt', 'delta_opt', ...
    'beta1_0', 'beta2_0', 'beta3_0', 'alpha0', ...
    'n', 'm_beta', 'delta_min', 'delta_max', ...
    'Np', 'phi', 'NumCurves');

fprintf('\n=== Results saved to Results/Poultry_Inverse_CA_tdep_INDEPENDENT_fix_results.mat ===\n');
fprintf('Estimated delta: %.3f (95%% CI: %.3f, %.3f)\n', delta_opt, prctile(Delta_params, 2.5), prctile(Delta_params, 97.5));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
function q = Phi_poultry(theta, tspan, X0, n, m_beta, phi, delta, Np, I_w, tw, a, b, ...
weight_incidence_1, weight_cumulative_1,weight_incidence_2, weight_cumulative_2,weight_incidence_3, weight_cumulative_3)
beta_params = reshape(theta(1:n*m_beta), m_beta, n)'; 
alpha = theta(n*m_beta+1:n*m_beta+6);

delta_vec = [delta; delta; delta];

alpha_matrix = build_alpha_matrix(alpha, n);

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[t, yp] = ode23s(@(t,y) farm_network(t, y, beta_params, n, phi, delta_vec, Np, ...
alpha_matrix, I_w, tw, a, b), tspan, X0, options);

S = yp(:,1:n);
I = yp(:,n+1:2*n);

q_incidence = zeros(length(t), n);
for k = 1:length(t)
Iw_t = interp1(tw, I_w, t(k), 'spline');
for i = 1:n
beta_i_t = beta(t(k), beta_params(i,:), a, b);
spatial = 0;
for j = 1:n
spatial = spatial + (I(k,j) / Np(j)) * alpha_matrix(i,j);
end
q_incidence(k,i) = beta_i_t * S(k,i) * Iw_t + S(k,i) * spatial;
end
end

q_cumulative = cumsum(q_incidence, 1);
q = [q_incidence(:,1)*weight_incidence_1; q_cumulative(:,1)*weight_cumulative_1; ...
q_incidence(:,2)*weight_incidence_2; q_cumulative(:,2)*weight_cumulative_2; ...
q_incidence(:,3)*weight_incidence_3; q_cumulative(:,3)*weight_cumulative_3];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = farm_network(t, y, beta_params, n, phi, delta, Np, alpha_matrix, I_w, tw, a, b)
S = y(1:n);
I = y(n+1:2*n);
C = y(2*n+1:3*n);
I_w_t = interp1(tw, I_w, t, 'spline');

dSdt = zeros(n,1);
dIdt = zeros(n,1);
dCdt = zeros(n,1);

for i = 1:n
beta_i_t = beta(t, beta_params(i,:), a, b);
spatial_infection = 0;
for j = 1:n
spatial_infection = spatial_infection + (I(j)/Np(j)) * alpha_matrix(i,j);
end
dSdt(i) = -beta_i_t*S(i)*I_w_t - S(i)*spatial_infection + delta(i)*C(i);
dIdt(i) = beta_i_t*S(i)*I_w_t + S(i)*spatial_infection - phi*I(i);
dCdt(i) = phi*I(i) - delta(i)*C(i);
end
dydt=[dSdt; dIdt; dCdt];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q = I_poultry(theta, tspan, X0, n, m_beta, phi, delta, Np, I_w, tw, a, b)
beta_params = reshape(theta(1:n*m_beta), m_beta, n)';
alpha = theta(n*m_beta+1:n*m_beta+6);

delta_vec =[delta; delta; delta];

alpha_matrix = build_alpha_matrix(alpha, n);

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[t, yp] = ode23s(@(t,y) farm_network(t, y, beta_params, n, phi, delta_vec, Np, ...
alpha_matrix, I_w, tw, a, b), tspan, X0, options);

S = yp(:,1:n);
I = yp(:,n+1:2*n);

q_incidence = zeros(length(t), n);
for k = 1:length(t)
Iw_t = interp1(tw, I_w, t(k), 'spline');
for i = 1:n
beta_i_t = beta(t(k), beta_params(i,:), a, b);
spatial = 0;
for j = 1:n
spatial = spatial+ (I(k,j) / Np(j)) * alpha_matrix(i,j);
end
q_incidence(k,i) =beta_i_t * S(k,i) * Iw_t + S(k,i) * spatial;
end
end

q_cumulative = cumsum(q_incidence, 1);
q = [q_incidence(:,1); q_cumulative(:,1); ...
q_incidence(:,2); q_cumulative(:,2); ...
q_incidence(:,3); q_cumulative(:,3)];
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function yy = beta(t, theta, a, b)
yy = 0;
m = length(theta);
for j = 1:m
yy = yy + theta(j)*leg(j-1,t,a,b);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = leg(j,t,a,b)
x = (2.*t - a - b)./(b - a);

if j == 0
P = 1;
elseif j == 1
P = x;
else
P1 = 1; P2 = x;
for k = 2:j
P3 = ((2*(k-1)+1).*x.*P2- (k-1).*P1)./k;
P1 = P2; P2 = P3;
end
P = P3;
end
end
