clc; clear; close all;
format long
warning('off','all')

% Needs Modification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
state = 'CA';
weight = 1;
t0 = 5*pi/6;
A1 = 117959600;
A2 = 215292600;
noise_lvl=10;
NumCurves = 100;
tdata = (1:38)';

% Loading
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wbird = readmatrix('wbird.xlsx');

DataI = wbird(:,2);
DataC = cumsum(DataI)/weight;
Data = [DataI;DataC];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=3;
a = min(tdata);
b = max(tdata);
N0 = (A1 + A2)/2 + (A2 - A1)/2 * sin((pi/6)*a + 4*pi/3-t0);
Ntdata = (A1 + A2)/2 + (A2 - A1)/2 * sin((pi/6)*tdata + 4*pi/3-t0);
% sigma_0 = (A1 + A2)*10^(-8)/2;
% sigma_0 = 1.5;
sigma_0 = .5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta0 = [sigma_0;3;1];

lb = [0;0;0];
ub = Inf * ones(m, 1);

 
curves_results = zeros(NumCurves, b-a+1);
incidence_results = zeros(NumCurves, b-a+1);
I_results = zeros(NumCurves, b-a+1);
Params_1 = zeros(NumCurves,m);

rng(0)
tau = repmat(abs(DataI(1:b-a+1)), [1, NumCurves]);
incidence_noisy = tau + normrnd(0,noise_lvl, size(tau));
negIdx = incidence_noisy < 0;
incidence_noisy(negIdx) = unifrnd(0, tau(negIdx));
curvesC = cumsum(incidence_noisy, 1); 
curves = [incidence_noisy;curvesC/weight];

parfor i = 1:NumCurves
try
ExpData = curves(:, i);
X0 = [N0 - ExpData(1); ExpData(1)];

options = optimoptions('lsqcurvefit', ...
'Algorithm', 'levenberg-marquardt', ...
'Display', 'iter', ...
'MaxIterations', 1000, ...
'TolX', 1e-10, ...
'TolFun', 1e-10, ...
'InitDamping', 1e7, ...
'UseParallel', true, ...
'SpecifyObjectiveGradient', false);

theta = lsqcurvefit(@(theta, tdata) Phi_wbird_ode(theta, A1, A2,t0, tdata, X0,weight), ...
theta0, tdata, ExpData, lb, ub, options);

[~, incidence_results(i, :), I_results(i, :)] = Phi_wbird_ode(theta, A1, A2,t0, tdata, X0,weight);
curves_results(i, :) = cumsum(incidence_results(i, :));
% theta(1) = 2*theta(1)/(A1 + A2);
Params_1(i, :) = theta;

catch
end
end

DataC_actual = DataC*weight;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
font_main     = 18;
font_axis     = 15;
font_legend   = 12;
font_date     = 14;
line_mean     = 1.5;
line_thin     = 1.0;
marker_size   = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure1 = figure;

if isnumeric(tdata)
    start_date = datetime(2022, 8, 1);
    tdata = start_date + calmonths(tdata - 1);
elseif iscell(tdata) || isstring(tdata)
    tdata = datetime(tdata);
end

t_start = tdata(1);
t_end   = tdata(end);
t_mid   = t_start + (t_end - t_start)/2;

subplot(1,3,1)
line1_I = plot(tdata, I_results', 'Color', [0, 0.6, 0, 0.1], 'LineWidth', line_thin, 'LineStyle', '-', 'HandleVisibility', 'off');
hold on
line2_I = plot(tdata, mean(I_results, 1), '-b', 'LineWidth', line_mean);
ylabel('Infected (Wild Birds)', 'Interpreter', 'latex', 'FontSize', font_main);
xlim([t_start t_end]);
xticks([t_start t_mid t_end]);
set(gca, 'XTickLabel', []);
text(t_start, min(ylim), datestr(t_start, 'mmm yyyy'), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
text(t_mid, min(ylim), datestr(t_mid, 'mmm yyyy'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
text(t_end, min(ylim), datestr(t_end, 'mmm yyyy'), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
legend([line1_I(1) line2_I(1)], {'Reconstructed Infected Curves', 'Mean of Reconstructed Infected'}, 'Interpreter', 'latex', 'Location', 'best', 'Color', 'none', 'FontSize', font_legend);
set(gca, 'FontSize', font_axis, 'TickLabelInterpreter', 'latex');

subplot(1,3,2)
line1_inc = plot(tdata, incidence_results, 'Color', [0.8, 0.4, 0, 0.1], 'LineWidth', line_thin, 'LineStyle', '-', 'HandleVisibility', 'off');
hold on
line2_inc = plot(tdata, mean(incidence_results), '-', 'Color', [0.8, 0.4, 0], 'LineWidth', line_mean, 'MarkerSize', marker_size);
line3_inc = plot(tdata, DataI, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', 'LineWidth', 1.2, 'MarkerSize', marker_size);
ylabel('Incidence (Wild Birds)', 'Interpreter', 'latex', 'FontSize', font_main);
xlim([t_start t_end]);
xticks([t_start t_mid t_end]);
set(gca, 'XTickLabel', []);
text(t_start, min(ylim), datestr(t_start, 'mmm yyyy'), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
text(t_mid, min(ylim), datestr(t_mid, 'mmm yyyy'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
text(t_end, min(ylim), datestr(t_end, 'mmm yyyy'), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
legend([line1_inc(1) line2_inc(1) line3_inc(1)], {'Reconstructed Incidence Curves', 'Mean of Reconstructed Incidence', 'Reported Incidence Data'}, 'Interpreter', 'latex', 'Location', 'best', 'Color', 'none', 'FontSize', font_legend);
set(gca, 'FontSize', font_axis, 'TickLabelInterpreter', 'latex');

subplot(1,3,3)
line1 = plot(tdata, curves_results, 'Color', [0.6, 0, 0, 0.1], 'LineWidth', line_thin, 'LineStyle', '-', 'HandleVisibility', 'off');
hold on
line2 = plot(tdata, mean(curves_results), '-r', 'LineWidth', line_mean, 'MarkerSize', marker_size);
line3 = plot(tdata, DataC, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', 'LineWidth', 1.2, 'MarkerSize', marker_size);
ylabel('Cumulative (Wild Birds)', 'Interpreter', 'latex', 'FontSize', font_main);
xlim([t_start t_end]);
xticks([t_start t_mid t_end]);
set(gca, 'XTickLabel', []);
text(t_start, min(ylim), datestr(t_start, 'mmm yyyy'), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
text(t_mid, min(ylim), datestr(t_mid, 'mmm yyyy'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
text(t_end, min(ylim), datestr(t_end, 'mmm yyyy'), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', font_date, 'Interpreter', 'latex');
legend([line1(1) line2(1) line3(1)], {'Reconstructed Cumulative Curves', 'Mean of Reconstructed Cumulative', 'Reported Cumulative Data'}, 'Interpreter', 'latex', 'Location', 'best', 'Color', 'none', 'FontSize', font_legend);
set(gca, 'FontSize', font_axis, 'TickLabelInterpreter', 'latex');

figure2 = figure;
param_names = {'\sigma', '\gamma_r', '\gamma_d'};
for j = 1:m
    subplot(1, m, j)
    h = histogram(Params_1(:, j), 'Normalization', 'probability');
    ylims = [0, max(h.Values)*1.1];
    ylim(ylims);
    hold on
    mean_val = mean(Params_1(:, j));
    ci_95 = prctile(Params_1(:, j), [2.5, 97.5]);
    plot([mean_val, mean_val], ylim, 'r-', 'LineWidth', 2)
    plot([ci_95(1), ci_95(1)], ylim, 'k--', 'LineWidth', 1.5)
    plot([ci_95(2), ci_95(2)], ylim, 'k--', 'LineWidth', 1.5)
    plot([theta0(j), theta0(j)], ylims, 'Color', [1 1 0], 'LineWidth', 2.2)
    title(sprintf('$%s = %.3f \\ (%.3f, %.3f)$', param_names{j}, mean_val, ci_95(1), ci_95(2)), 'Interpreter', 'latex', 'FontSize', font_main);
    set(gca, 'FontSize', font_axis, 'TickLabelInterpreter', 'latex');
end

exportFig(figure1, 'panel', 15, 4, state);
exportFig(figure2, 'params', 15, 4, state);


save(sprintf('wildbird_%s.mat', state), ...
    'tdata', ...
    'curves_results', ...
    'DataC', ...
    'incidence_results', ...
    'DataI', ...
    'I_results', ...
    'Params_1', ...
    'theta0', ...
    'm', ...
    'Ntdata');

mean_I = mean(I_results, 1);
T = table(tdata(:), mean_I(:), ...
    'VariableNames', {'Time', 'Mean_Infected'});
writetable(T, sprintf('wildbird_Iw_%s.xlsx', state));

T_all = table(tdata(:));
T_all.Properties.VariableNames{1} = 'Time';
for i = 1:NumCurves
    T_all.(sprintf('I_iter_%d', i)) = I_results(i, :)';
end
writetable(T_all, sprintf('wildbird_Iw_%s_all.xlsx', state));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = wbird_ode(t, y, sigma, gamma_wr, gamma_wd, A1, A2, t0)
Nprime = (A2 - A1)/2 * (pi/6) * cos((pi/6)*t  + 4*pi/3-t0);
N_a = (A1 + A2)/2;
dydt = zeros(2,1);
dydt(1) = -sigma*y(1)*y(2)/N_a+ gamma_wr*y(2) + Nprime;   
dydt(2) =  sigma*y(1)*y(2)/N_a- gamma_wr*y(2) - gamma_wd*y(2); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [q, incidence, I] = Phi_wbird_ode(theta,A1,A2,t0, tspan,X0,weight)
sigma = theta(1);
gamma_wr = theta(2);
gamma_wd = theta(3);
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[~, Y] = ode45(@(t, y) wbird_ode(t, y, sigma, gamma_wr, gamma_wd, A1, A2,t0), tspan, X0, options);
S = Y(:,1);
I = Y(:,2);
N_a = (A1 + A2)/2;
incidence = sigma.*S.*I./N_a;
q1 = cumsum(incidence);
q = [incidence;q1/weight];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function exportFig(fig, name, width, height, state)
set(fig,'Units','inches','Position',[0 0 width height]);
set(fig,'PaperUnits','inches','PaperPosition',[0 0 width height],'PaperSize',[width height]);
exportgraphics(fig, sprintf('wildbird_%s_%s.pdf', name, state), 'ContentType','vector');
end


