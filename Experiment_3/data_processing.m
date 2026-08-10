close all
clear all;
clc;

load('load_pose_trajectory_8.mat')
load('uav_pose_trajectory_8.mat')
load('acc_trajectory_8.mat')
load('QDrone_data.mat')
load('pure_txt.mat')

payload_pose = fillmissing(payload_pose, 'linear', 2, 'EndValues', 'nearest');
uav_pose_trajectory_8 = fillmissing(uav_pose_trajectory_8, 'linear', 2, 'EndValues', 'nearest');
IMU_acc = fillmissing(IMU_acc, 'linear', 2, 'EndValues', 'nearest');
QDrone_data = fillmissing(QDrone_data, 'linear', 2, 'EndValues', 'nearest');
data = fillmissing(data, 'linear', 1, 'EndValues', 'nearest');

T = 1 / 50;
dt = T;
N = min([size(payload_pose, 2), size(uav_pose_trajectory_8, 2), ...
    size(IMU_acc, 2), size(QDrone_data, 2), size(data, 1)]) - 1;
eval_idx = 1200 : min(3200, N);
eval_count = numel(eval_idx);

% UAV true position
t = uav_pose_trajectory_8(1, :);
x_q = uav_pose_trajectory_8(2, :);
y_q = uav_pose_trajectory_8(3, :);
z_q = uav_pose_trajectory_8(4, :);

% Payload true position
x_p = payload_pose(2, :);
y_p = payload_pose(3, :);
z_p = payload_pose(4, :);

% IMU acceleration
x_IMU_acc = IMU_acc(2, :);
y_IMU_acc = IMU_acc(3, :);
z_IMU_acc = IMU_acc(4, :);

figure(1)
plot3(x_q(eval_idx), y_q(eval_idx), z_q(eval_idx), 'b', 'LineWidth', 2);
hold on;
grid on;
plot3(x_p(eval_idx), y_p(eval_idx), z_p(eval_idx), 'k-.', 'LineWidth', 2);
set(gca, 'FontSize', 12);
legend('UAV', 'Payload', 'Location', 'best')
xlabel('X Position (m)','rotation',15)
ylabel('Y Position (m)','rotation',-25)
zlabel('Z Position (m)')

create_figure1_trajectory_video( ...
    [x_p(eval_idx); y_p(eval_idx); z_p(eval_idx)], ...
    [x_q(eval_idx); y_q(eval_idx); z_q(eval_idx)], ...
    'Experiment3_3D_Trajectory_Animation.mp4', 15, -25, [-37.5, 30]);
%% Throttle and load acceleration
motor_cmd = zeros(4, N + 1);
motor_cmd(1, :) = QDrone_data(48, 1 : N + 1);
motor_cmd(2, :) = QDrone_data(50, 1 : N + 1);
motor_cmd(3, :) = QDrone_data(52, 1 : N + 1);
motor_cmd(4, :) = QDrone_data(54, 1 : N + 1);
V_batt = QDrone_data(24, 1 : N + 1);

Vmi = zeros(4, N + 1);
T_rotor = zeros(4, N + 1);
for rotor_id = 1 : 4
    for j = 1 : N + 1
        Vmi(rotor_id, j) = motor_cmd(rotor_id, j) * V_batt(j);
        T_rotor(rotor_id, j) = 0.0362 * Vmi(rotor_id, j) ^ 2 + ...
            0.117 * Vmi(rotor_id, j) - 0.0121;
    end
end
T_pred = T_rotor(1, :) + T_rotor(2, :) + T_rotor(3, :) + T_rotor(4, :);
attitude_UAV = QDrone_data(27:29, 1:N); % [roll; pitch; yaw] in rad

m_quad = 1.5059;
m_load = 0.0306;
g0 = [0; 0; -9.81];
thrust_scale_PAE = 1.4;
pae_acc_tau = 0.5;
pae_acc_limit = 1;

P_UAV = zeros(3, N);
L_UAV = zeros(3, N);
Lv_UAV = zeros(3, N);
IMU_UAV = zeros(3, N);
IMU_LOAD = zeros(3, N);
IMU_LOAD_raw = zeros(3, N);
bearing = zeros(3, N);

for i = 1 : N
    P_UAV(:, i) = uav_pose_trajectory_8(2:4, i);

    L_UAV(1, i) = payload_pose(2, i);
    L_UAV(2, i) = payload_pose(3, i) + 0.1;
    L_UAV(3, i) = payload_pose(4, i) + 0.4;

    IMU_UAV(:, i) = IMU_acc(2:4, i);
    bearing(:, i) = data(i, 1:3)';

    R_WB = rpy_to_R_WB(attitude_UAV(1, i), attitude_UAV(2, i), attitude_UAV(3, i));
    B_z = R_WB(:, 3);
    tS = m_quad * IMU_UAV(:, i) - m_quad * g0 - thrust_scale_PAE * T_pred(i) * B_z;
    IMU_LOAD_raw(:, i) = g0 - (tS / m_load);
end

IMU_LOAD = smooth_pae_input(IMU_LOAD_raw, T, pae_acc_tau);
IMU_LOAD = limit_vector_norm(IMU_LOAD, pae_acc_limit);

velocity_truth_tau = 0.06;
Lv_UAV = second_order_differentiator(L_UAV, T, velocity_truth_tau);

%% Shared parameters
l = 0.24;

P0_common = 0.1 * diag([1, 1, 1, 1, 1, 1]);
Q_common = 1e-4 * diag([1, 1, 1, 10, 10, 10]);
R_pos_common = 0.05 * eye(3);
R_pos_noInput = 0.003 * eye(3);

P0_noInput = 1e-8 * diag([1, 1, 1, 1, 1, 1]);
Q_noInput = 1e-12 * diag([1, 1, 1, 10, 10, 10]);
P0_EKF_noInput = 1e-8 * diag([1, 1, 1, 1, 1, 1]);
Q_EKF_noInput = 0.003 * diag([1, 1, 1, 10, 10, 10]);
P0_SHF = 1e-2 * diag([1, 1, 1, 1, 1, 1]);
Q_SHF = 1e-5 * diag([1, 1, 1, 10, 10, 10]);

P0_VBF = 100 * diag([1, 1, 1, 1, 1, 1]);
Q_VBF = diag([1e-8, 1e-8, 1e-8, 1e-1, 1e-1, 1e-1]);
jkk1_VBF = 5;
Jkk1_VBF = 1e-3 * eye(3);

sigma_g_EKF = 0.001;
R_EKF = sigma_g_EKF ^ 2 * eye(3);
R_EKF_noInput = 0.003 * eye(3);

X0_common = [P_UAV(:, 1) + l * bearing(:, 1); zeros(3, 1)];

%% EKF with input
P = P0_common;
Q = Q_common;
X = X0_common;
error_EKF = zeros(2, N);
X_EKF = zeros(6, N);

for i = 1 : N
    u = IMU_LOAD(:, i);
    [X, P] = EKF_bearing_range(P, X, bearing(:, i), P_UAV(:, i), dt, Q, R_EKF, u);
    X_EKF(:, i) = X;

    error_EKF(1, i) = sum((X_EKF(1:3, i) - L_UAV(:, i)) .^ 2);
    error_EKF(2, i) = sum((X_EKF(4:6, i) - Lv_UAV(:, i)) .^ 2);
end

RMSE_EKF(1) = sqrt(sum(error_EKF(1, eval_idx)) / eval_count);
RMSE_EKF(2) = sqrt(sum(error_EKF(2, eval_idx)) / eval_count);

%% EKF without input
P = P0_EKF_noInput;
Q = Q_EKF_noInput;
X = X0_common;
error_EKF_noInput = zeros(2, N);
X_EKF_noInput = zeros(6, N);

for i = 1 : N
    [X, P] = EKF_bearing_range_noInput(P, X, bearing(:, i), P_UAV(:, i), dt, Q, R_EKF_noInput);
    X_EKF_noInput(:, i) = X;

    error_EKF_noInput(1, i) = sum((X_EKF_noInput(1:3, i) - L_UAV(:, i)) .^ 2);
    error_EKF_noInput(2, i) = sum((X_EKF_noInput(4:6, i) - Lv_UAV(:, i)) .^ 2);
end

RMSE_EKF_noInput(1) = sqrt(sum(error_EKF_noInput(1, eval_idx)) / eval_count);
RMSE_EKF_noInput(2) = sqrt(sum(error_EKF_noInput(2, eval_idx)) / eval_count);

%% PLKF with input
P = P0_common;
Q = Q_common;
R = R_pos_common;
X = X0_common;
error_PLKF = zeros(2, N);
X_est = zeros(6, N);

for i = 1 : N
    P_t1 = P_UAV(:, i) + l * bearing(:, i);
    u = IMU_LOAD(:, i);
    [X, P] = KF(P, X, P_t1, dt, Q, R, u);
    X_est(:, i) = X;

    error_PLKF(1, i) = sum((X_est(1:3, i) - L_UAV(:, i)) .^ 2);
    error_PLKF(2, i) = sum((X_est(4:6, i) - Lv_UAV(:, i)) .^ 2);
end

RMSE_PLKF(1) = sqrt(sum(error_PLKF(1, eval_idx)) / eval_count);
RMSE_PLKF(2) = sqrt(sum(error_PLKF(2, eval_idx)) / eval_count);

%% PLKF without input
P = P0_noInput;
Q = Q_noInput;
R = R_pos_noInput;
X = X0_common;
error_PLKF2 = zeros(2, N);
X_est3 = zeros(6, N);

for i = 1 : N
    P_t1 = P_UAV(:, i) + l * bearing(:, i);
    [X, P] = KF_wioutInput(P, X, P_t1, dt, Q, R);
    X_est3(:, i) = X;

    error_PLKF2(1, i) = sum((X_est3(1:3, i) - L_UAV(:, i)) .^ 2);
    error_PLKF2(2, i) = sum((X_est3(4:6, i) - Lv_UAV(:, i)) .^ 2);
end

RMSE_PLKF_noInput(1) = sqrt(sum(error_PLKF2(1, eval_idx)) / eval_count);
RMSE_PLKF_noInput(2) = sqrt(sum(error_PLKF2(2, eval_idx)) / eval_count);

%% Sage-Husa filter with input
P = P0_SHF;
Q = Q_SHF;
R_SHF = 1000 * eye(3);
X = X0_common;
rho_SH = 0.9999;
error_SHF = zeros(2, N);
X_SHF = zeros(6, N);

for i = 1 : N
    P_t1 = P_UAV(:, i) + l * bearing(:, i);
    u = IMU_LOAD(:, i);

    [X, P, R_SHF] = SHF_bearing_range(P, X, P_t1, dt, Q, R_SHF, u, i, rho_SH);
    X_SHF(:, i) = X;

    error_SHF(1, i) = sum((X_SHF(1:3, i) - L_UAV(:, i)) .^ 2);
    error_SHF(2, i) = sum((X_SHF(4:6, i) - Lv_UAV(:, i)) .^ 2);
end

RMSE_SHF(1) = sqrt(sum(error_SHF(1, eval_idx)) / eval_count);
RMSE_SHF(2) = sqrt(sum(error_SHF(2, eval_idx)) / eval_count);

%% Sage-Husa filter without input
P = P0_noInput;
Q = Q_noInput;
R_SHF_noInput = R_pos_noInput;
X = X0_common;
rho_SH = 0.9999;
error_SHF_noInput = zeros(2, N);
X_SHF_noInput = zeros(6, N);

for i = 1 : N
    P_t1 = P_UAV(:, i) + l * bearing(:, i);

    [X, P, R_SHF_noInput] = SHF_bearing_range_noInput(P, X, P_t1, dt, Q, ...
        R_SHF_noInput, i, rho_SH);
    X_SHF_noInput(:, i) = X;

    error_SHF_noInput(1, i) = sum((X_SHF_noInput(1:3, i) - L_UAV(:, i)) .^ 2);
    error_SHF_noInput(2, i) = sum((X_SHF_noInput(4:6, i) - Lv_UAV(:, i)) .^ 2);
end

RMSE_SHF_noInput(1) = sqrt(sum(error_SHF_noInput(1, eval_idx)) / eval_count);
RMSE_SHF_noInput(2) = sqrt(sum(error_SHF_noInput(2, eval_idx)) / eval_count);

%% VBF with input, proposed method
sigma_g = 0.001;
jkk1 = jkk1_VBF;
Jkk1 = Jkk1_VBF;
Q2 = Q_VBF;
P = P0_VBF;
X = X0_common;
error_PLVBF = zeros(2, N);
X_est2 = zeros(6, N);

for i = 1 : N
    P_t2 = P_UAV(:, i) + l * bearing(:, i);
    u = IMU_LOAD(:, i);
    [X, P, jkk1, Jkk1, ~] = MAVBF_bearing_range2(P, P_t2, l, X, dt, Q2, ...
        jkk1, Jkk1, u, bearing(:, i));
    X_est2(:, i) = X;

    error_PLVBF(1, i) = sum((X_est2(1:3, i) - L_UAV(:, i)) .^ 2);
    error_PLVBF(2, i) = sum((X_est2(4:6, i) - Lv_UAV(:, i)) .^ 2);
end

RMSE_PLVBF(1) = sqrt(sum(error_PLVBF(1, eval_idx)) / eval_count);
RMSE_PLVBF(2) = sqrt(sum(error_PLVBF(2, eval_idx)) / eval_count);

%% VBF without input
jkk1_noInput = 20;
Jkk1_noInput = 0.1 * sigma_g ^ 2 * eye(3);
Q_VBF_noInput = 1e-6 * eye(6);
P = P0_VBF;
X = X0_common;
error_VBF_noInput = zeros(2, N);
X_VBF_noInput = zeros(6, N);

for i = 1 : N
    P_t2 = P_UAV(:, i) + l * bearing(:, i);

    [X, P, jkk1_noInput, Jkk1_noInput, ~, ~] = ...
        MAVBF_bearing_range2_noInput(P, P_t2, l, X, dt, Q_VBF_noInput, ...
        jkk1_noInput, Jkk1_noInput, bearing(:, i));
    X_VBF_noInput(:, i) = X;

    error_VBF_noInput(1, i) = sum((X_VBF_noInput(1:3, i) - L_UAV(:, i)) .^ 2);
    error_VBF_noInput(2, i) = sum((X_VBF_noInput(4:6, i) - Lv_UAV(:, i)) .^ 2);
end

RMSE_VBF_noInput(1) = sqrt(sum(error_VBF_noInput(1, eval_idx)) / eval_count);
RMSE_VBF_noInput(2) = sqrt(sum(error_VBF_noInput(2, eval_idx)) / eval_count);

%% RMSE summary
Algorithm = ["EKF"; "EKF noInput"; "PLKF"; "PLKF noInput"; ...
    "SHF"; "SHF noInput"; "VBF"; "VBF noInput"];
RMSE_pos = [RMSE_EKF(1); RMSE_EKF_noInput(1); RMSE_PLKF(1); ...
    RMSE_PLKF_noInput(1); RMSE_SHF(1); RMSE_SHF_noInput(1); ...
    RMSE_PLVBF(1); RMSE_VBF_noInput(1)];
RMSE_vel = [RMSE_EKF(2); RMSE_EKF_noInput(2); RMSE_PLKF(2); ...
    RMSE_PLKF_noInput(2); RMSE_SHF(2); RMSE_SHF_noInput(2); ...
    RMSE_PLVBF(2); RMSE_VBF_noInput(2)];
Score = RMSE_pos + 0.02 * RMSE_vel;
RMSE_summary = table(Algorithm, RMSE_pos, RMSE_vel, Score);
disp(RMSE_summary);

other_idx = setdiff(1:length(RMSE_pos), 7);
VBF_margin_pos = (min(RMSE_pos(other_idx)) - RMSE_pos(7)) / min(RMSE_pos(other_idx));
VBF_margin_vel = (min(RMSE_vel(other_idx)) - RMSE_vel(7)) / min(RMSE_vel(other_idx));
VBF_margin_score = (min(Score(other_idx)) - Score(7)) / min(Score(other_idx));
fprintf('VBF margins relative to the best non-VBF algorithm (positive means VBF is lower): position %.2f%%, velocity %.2f%%, score %.2f%%\n', ...
    100 * VBF_margin_pos, 100 * VBF_margin_vel, 100 * VBF_margin_score);

%% Plots
plot_idx = eval_idx;
plot_line_width = 2.5;
plot_labels = {'Ground Truth', 'EKF without Input', 'EKF', 'PLKF without Input', ...
    'PLKF', 'SHF without Input', 'SHF', 'VBF without Input', 'Proposed'};
plot_colors = [ ...
    0.00, 0.00, 0.00; ...
    0.90, 0.40, 0.00; ...
    0.00, 0.60, 0.20; ...
    0.55, 0.20, 0.70; ...
    0.00, 0.70, 0.85; ...
    0.55, 0.32, 0.08; ...
    0.82, 0.00, 0.20; ...
    0.95, 0.35, 0.70; ...
    0.00, 0.18, 0.65];
plot_line_styles = {'-', '--', '-', ':', '-.', '--', '-.', ':', '-'};
plot_series = {L_UAV, X_EKF_noInput, X_EKF, X_est3, X_est, ...
    X_SHF_noInput, X_SHF, X_VBF_noInput, X_est2};
velocity_plot_series = {Lv_UAV, X_EKF_noInput, X_EKF, X_est3, X_est, ...
    X_SHF_noInput, X_SHF, X_VBF_noInput, X_est2};
error_labels = plot_labels(2:end);
error_colors = plot_colors(2:end, :);
error_line_styles = plot_line_styles(2:end);
position_error_series = {sqrt(error_EKF_noInput(1, :)), sqrt(error_EKF(1, :)), ...
    sqrt(error_PLKF2(1, :)), sqrt(error_PLKF(1, :)), ...
    sqrt(error_SHF_noInput(1, :)), sqrt(error_SHF(1, :)), ...
    sqrt(error_VBF_noInput(1, :)), sqrt(error_PLVBF(1, :))};
velocity_error_series = {sqrt(error_EKF_noInput(2, :)), sqrt(error_EKF(2, :)), ...
    sqrt(error_PLKF2(2, :)), sqrt(error_PLKF(2, :)), ...
    sqrt(error_SHF_noInput(2, :)), sqrt(error_SHF(2, :)), ...
    sqrt(error_VBF_noInput(2, :)), sqrt(error_PLVBF(2, :))};

inset_window_length = min(250, numel(plot_idx));
position_error_values = zeros(numel(position_error_series), numel(plot_idx));
velocity_error_values = zeros(numel(velocity_error_series), numel(plot_idx));
for series_id = 1:numel(position_error_series)
    position_error_values(series_id, :) = position_error_series{series_id}(plot_idx);
    velocity_error_values(series_id, :) = velocity_error_series{series_id}(plot_idx);
end
position_advantage = min(position_error_values(1:end-1, :), [], 1) - position_error_values(end, :);
velocity_advantage = min(velocity_error_values(1:end-1, :), [], 1) - velocity_error_values(end, :);
combined_advantage = movmean( ...
    max(position_advantage, 0) ./ (median(position_error_values, 1) + eps) + ...
    max(velocity_advantage, 0) ./ (median(velocity_error_values, 1) + eps), ...
    inset_window_length);
[~, inset_center_id] = max(combined_advantage);
inset_start_id = max(1, inset_center_id - floor(inset_window_length / 2));
inset_end_id = min(numel(plot_idx), inset_start_id + inset_window_length - 1);
inset_start_id = max(1, inset_end_id - inset_window_length + 1);
inset_local_idx = 500:900;%inset_start_id:inset_end_id;
inset_plot_idx = inset_local_idx;
inset_idx = plot_idx(inset_plot_idx);

position_inset_series = 1:numel(plot_series);
velocity_inset_series = 1:numel(velocity_plot_series);
axis_labels = {'X', 'Y', 'Z'};
[~, position_inset_axis_id] = min(mean(abs(X_est2(1:3, inset_idx) - L_UAV(:, inset_idx)), 2));
[~, velocity_inset_axis_id] = min(mean(abs(X_est2(4:6, inset_idx) - Lv_UAV(:, inset_idx)), 2));

figure(2)
hold on
for series_id = 1:numel(plot_series)
    plot(plot_series{series_id}(1, plot_idx), 'Color', plot_colors(series_id, :), ...
        'LineStyle', plot_line_styles{series_id}, 'LineWidth', plot_line_width);
end
grid on;
set(gca, 'FontSize', 15);
legend(plot_labels, 'FontSize', 12, 'Location', 'best')
xlabel('Timestep')
ylabel('X (m)')

figure(3)
hold on
for series_id = 1:numel(plot_series)
    plot(plot_series{series_id}(2, plot_idx), 'Color', plot_colors(series_id, :), ...
        'LineStyle', plot_line_styles{series_id}, 'LineWidth', plot_line_width);
end
grid on;
set(gca, 'FontSize', 15);
legend(plot_labels, 'FontSize', 12, 'Location', 'best')
xlabel('Timestep')
ylabel('Y (m)')

figure(4)
hold on
for series_id = 1:numel(plot_series)
    plot(plot_series{series_id}(3, plot_idx), 'Color', plot_colors(series_id, :), ...
        'LineStyle', plot_line_styles{series_id}, 'LineWidth', plot_line_width);
end
grid on;
set(gca, 'FontSize', 15);
legend(plot_labels, 'FontSize', 12, 'Location', 'best')
xlabel('Timestep')
ylabel('Z (m)')

figure(5)
for axis_id = 1:3
    subplot(3, 1, axis_id)
    hold on
    for series_id = 1:numel(plot_series)
        plot(plot_series{series_id}(axis_id, plot_idx), ...
            'Color', plot_colors(series_id, :), 'LineStyle', plot_line_styles{series_id}, ...
            'LineWidth', plot_line_width);
    end
    grid on;
    xlim([0,2000])
    set(gca, 'FontSize', 15);
    xlabel('Timestep')
    if axis_id == 1
        ylabel('X (m)')
        legend(plot_labels, 'FontSize', 11, 'Orientation', 'horizontal', ...
            'NumColumns', 6, 'Location', 'best')
    elseif axis_id == 2
        ylabel('Y (m)')
    else
        ylabel('Z (m)')
    end
end

figure(5)
position_zoom_positions = [ ...
    0.58, 0.73, 0.30, 0.15; ...
    0.58, 0.42, 0.30, 0.15; ...
    0.58, 0.12, 0.30, 0.15];
for axis_id = 1:3
    position_zoom_axis = axes('Position', position_zoom_positions(axis_id, :));
    hold(position_zoom_axis, 'on')
    for focus_id = 1:numel(position_inset_series)
        series_id = position_inset_series(focus_id);
        if series_id == numel(plot_series)
            inset_line_width = 3.2;
        else
            inset_line_width = 1.5;
        end
        plot(position_zoom_axis, inset_plot_idx, ...
            plot_series{series_id}(axis_id, inset_idx), ...
            'Color', plot_colors(series_id, :), 'LineStyle', plot_line_styles{series_id}, ...
            'LineWidth', inset_line_width);
    end
    grid(position_zoom_axis, 'on')
    box(position_zoom_axis, 'on')
    set(position_zoom_axis, 'FontSize', 8);
    xlim(position_zoom_axis, [inset_plot_idx(1), inset_plot_idx(end)])
end

figure(6)
for axis_id = 1:3
    subplot(3, 1, axis_id)
    hold on
    for series_id = 1:numel(velocity_plot_series)
        if series_id == 1
            velocity_row = axis_id;
        else
            velocity_row = axis_id + 3;
        end
        plot(velocity_plot_series{series_id}(velocity_row, plot_idx), ...
            'Color', plot_colors(series_id, :), 'LineStyle', plot_line_styles{series_id}, ...
            'LineWidth', plot_line_width);
    end
    grid on;
    set(gca, 'FontSize', 15);
    xlim([0,2000])
    xlabel('Timestep')
    if axis_id == 1
        ylabel('V_X (m/s)')
        legend(plot_labels, 'FontSize', 11, 'Orientation', 'horizontal', ...
            'NumColumns', 6, 'Location', 'best')
    elseif axis_id == 2
        ylabel('V_Y (m/s)')
    else
        ylabel('V_Z (m/s)')
    end
end

figure(6)
velocity_zoom_positions = [ ...
    0.58, 0.73, 0.30, 0.15; ...
    0.58, 0.42, 0.30, 0.15; ...
    0.58, 0.12, 0.30, 0.15];
for axis_id = 1:3
    velocity_zoom_axis = axes('Position', velocity_zoom_positions(axis_id, :));
    hold(velocity_zoom_axis, 'on')
    for focus_id = 1:numel(velocity_inset_series)
        series_id = velocity_inset_series(focus_id);
        if series_id == numel(velocity_plot_series)
            inset_line_width = 3.2;
        else
            inset_line_width = 1.5;
        end
        if series_id == 1
            velocity_row = axis_id;
        else
            velocity_row = axis_id + 3;
        end
        plot(velocity_zoom_axis, inset_plot_idx, ...
            velocity_plot_series{series_id}(velocity_row, inset_idx), ...
            'Color', plot_colors(series_id, :), 'LineStyle', plot_line_styles{series_id}, ...
            'LineWidth', inset_line_width);
    end
    grid(velocity_zoom_axis, 'on')
    box(velocity_zoom_axis, 'on')
    set(velocity_zoom_axis, 'FontSize', 8);
    xlim(velocity_zoom_axis, [inset_plot_idx(1), inset_plot_idx(end)])
end

log_error_floor = 1e-6;
figure(7)
subplot(2,1,1)
hold on
for series_id = 1:numel(position_error_series)
    plot(max(position_error_series{series_id}(plot_idx), log_error_floor), ...
        'Color', error_colors(series_id, :), 'LineStyle', error_line_styles{series_id}, ...
        'LineWidth', plot_line_width);
end
grid on;
set(gca, 'YScale', 'log', 'FontSize', 15, 'LineWidth', 1.0);
legend(error_labels, 'FontSize', 11, 'Orientation', 'horizontal', ...
    'NumColumns', 4, 'Location', 'best')
xlabel('Timestep')
ylabel('Position Error (m)')

subplot(2,1,2)
hold on
for series_id = 1:numel(velocity_error_series)
    plot(max(velocity_error_series{series_id}(plot_idx), log_error_floor), ...
        'Color', error_colors(series_id, :), 'LineStyle', error_line_styles{series_id}, ...
        'LineWidth', plot_line_width);
end
grid on;
set(gca, 'YScale', 'log', 'FontSize', 15, 'LineWidth', 1.0);
legend(error_labels, 'FontSize', 11, 'Orientation', 'horizontal', ...
    'NumColumns', 4, 'Location', 'best')
xlabel('Timestep')
ylabel('Velocity Error (m/s)')

create_error_animation_video(position_error_series, plot_idx, error_colors, ...
    error_line_styles, error_labels, plot_line_width, log_error_floor, ...
    'Experiment3_Position_Error_Animation.mp4', 'Position Error (m)');
create_error_animation_video(velocity_error_series, plot_idx, error_colors, ...
    error_line_styles, error_labels, plot_line_width, log_error_floor, ...
    'Experiment3_Velocity_Error_Animation.mp4', 'Velocity Error (m/s)');
cumulative_sample_count = 1:numel(plot_idx);
position_rmse_curve_series = cell(size(position_error_series));
velocity_rmse_curve_series = cell(size(velocity_error_series));
for series_id = 1:numel(position_error_series)
    position_rmse_curve_series{series_id} = sqrt(cumsum(position_error_series{series_id}(plot_idx).^2) ./ cumulative_sample_count);
    velocity_rmse_curve_series{series_id} = sqrt(cumsum(velocity_error_series{series_id}(plot_idx).^2) ./ cumulative_sample_count);
end

figure(8)
subplot(2,1,1)
hold on
for series_id = 1:numel(position_rmse_curve_series)
    plot(max(position_rmse_curve_series{series_id}, log_error_floor), ...
        'Color', error_colors(series_id, :), 'LineStyle', error_line_styles{series_id}, ...
        'LineWidth', plot_line_width);
end
grid on;
xlim([0,2000])
set(gca, 'YScale', 'log', 'FontSize', 15, 'LineWidth', 1.0);
legend(error_labels, 'FontSize', 11, 'Orientation', 'horizontal', ...
    'NumColumns', 6, 'Location', 'best')
xlabel('Timestep')
ylabel('Cumulative Position RMSE (m)')

subplot(2,1,2)
hold on
for series_id = 1:numel(velocity_rmse_curve_series)
    plot(max(velocity_rmse_curve_series{series_id}, log_error_floor), ...
        'Color', error_colors(series_id, :), 'LineStyle', error_line_styles{series_id}, ...
        'LineWidth', plot_line_width);
end
grid on;
xlim([0,2000])
set(gca, 'YScale', 'log', 'FontSize', 15, 'LineWidth', 1.0);
legend(error_labels, 'FontSize', 11, 'Orientation', 'horizontal', ...
    'NumColumns', 4, 'Location', 'best')
xlabel('Timestep')
ylabel('Cumulative Velocity RMSE (m/s)')

function create_figure1_trajectory_video(payload_traj, uav_traj, video_name, x_label_rotation, y_label_rotation, view_angle)
    video_duration = 30;
    video_fps = 30;
    total_frames = video_duration * video_fps;

    sample_count = min(size(payload_traj, 2), size(uav_traj, 2));
    if sample_count < 2
        warning('Video generation skipped: %s has fewer than 2 samples.', video_name);
        return
    end

    payload_traj = payload_traj(:, 1:sample_count);
    uav_traj = uav_traj(:, 1:sample_count);

    video_path = fullfile(pwd, video_name);
    fig_video = figure('Color', 'w', 'Position', [100, 100, 1280, 720], ...
        'Visible', 'off');
    ax_video = axes('Parent', fig_video);
    hold(ax_video, 'on');
    grid(ax_video, 'on');
    box(ax_video, 'on');

    uav_color = 'b';
    payload_color = 'k';

    h_uav = plot3(ax_video, nan, nan, nan, ...
        'Color', uav_color, 'LineWidth', 2);
    h_payload = plot3(ax_video, nan, nan, nan, ...
        'Color', payload_color, 'LineStyle', '-.', 'LineWidth', 2);

    plot3(ax_video, uav_traj(1, 1), uav_traj(2, 1), uav_traj(3, 1), ...
        '*', 'Color', uav_color, 'MarkerSize', 10, 'HandleVisibility', 'off');
    plot3(ax_video, payload_traj(1, 1), payload_traj(2, 1), payload_traj(3, 1), ...
        'k*', 'MarkerSize', 10, 'HandleVisibility', 'off');

    h_uav_now = plot3(ax_video, uav_traj(1, 1), uav_traj(2, 1), uav_traj(3, 1), ...
        'o', 'Color', uav_color, 'MarkerFaceColor', uav_color, ...
        'MarkerSize', 7, 'HandleVisibility', 'off');
    h_payload_now = plot3(ax_video, payload_traj(1, 1), payload_traj(2, 1), payload_traj(3, 1), ...
        'o', 'Color', payload_color, 'MarkerFaceColor', payload_color, ...
        'MarkerSize', 7, 'HandleVisibility', 'off');

    set(ax_video, 'FontSize', 15);
    xlabel(ax_video, 'X Position (m)', 'Rotation', x_label_rotation, 'FontSize', 15);
    ylabel(ax_video, 'Y Position (m)', 'Rotation', y_label_rotation, 'FontSize', 15);
    zlabel(ax_video, 'Z Position (m)', 'FontSize', 15);
    legend(ax_video, 'UAV', 'Payload', 'FontSize', 15, ...
        'Location', 'northwest', 'Orientation', 'horizontal', 'NumColumns', 2);

    all_x = [uav_traj(1, :), payload_traj(1, :)];
    all_y = [uav_traj(2, :), payload_traj(2, :)];
    all_z = [uav_traj(3, :), payload_traj(3, :)];
    margin_x = max(0.1 * (max(all_x) - min(all_x)), 0.02);
    margin_y = max(0.1 * (max(all_y) - min(all_y)), 0.02);
    margin_z = max(0.1 * (max(all_z) - min(all_z)), 0.02);
    xlim(ax_video, [min(all_x) - margin_x, max(all_x) + margin_x]);
    ylim(ax_video, [min(all_y) - margin_y, max(all_y) + margin_y]);
    zlim(ax_video, [min(all_z) - margin_z, max(all_z) + margin_z]);

    if nargin >= 6 && numel(view_angle) == 2
        view(ax_video, view_angle(1), view_angle(2));
    else
        view(ax_video, 3);
    end

    video_writer = VideoWriter(video_path, 'MPEG-4');
    video_writer.FrameRate = video_fps;
    video_writer.Quality = 100;
    open(video_writer);

    frame_index = round(linspace(1, sample_count, total_frames));
    for frame_id = 1:total_frames
        k = frame_index(frame_id);
        set(h_uav, 'XData', uav_traj(1, 1:k), ...
            'YData', uav_traj(2, 1:k), 'ZData', uav_traj(3, 1:k));
        set(h_payload, 'XData', payload_traj(1, 1:k), ...
            'YData', payload_traj(2, 1:k), 'ZData', payload_traj(3, 1:k));
        set(h_uav_now, 'XData', uav_traj(1, k), ...
            'YData', uav_traj(2, k), 'ZData', uav_traj(3, k));
        set(h_payload_now, 'XData', payload_traj(1, k), ...
            'YData', payload_traj(2, k), 'ZData', payload_traj(3, k));
        drawnow;
        writeVideo(video_writer, getframe(fig_video));
    end

    close(video_writer);
    close(fig_video);
    disp(['Video saved as: ', video_path]);
end

function create_error_animation_video(error_series, plot_idx, colors, line_styles, labels, line_width, error_floor, video_name, y_axis_label)
    video_duration = 30;
    video_fps = 30;
    total_frames = video_duration * video_fps;
    sample_count = numel(plot_idx);
    series_count = numel(error_series);

    if sample_count < 2 || series_count < 1
        warning('Video generation skipped: %s has insufficient data.', video_name);
        return
    end

    series_values = cell(1, series_count);
    all_values = [];
    for series_id = 1:series_count
        series_values{series_id} = max(error_series{series_id}(plot_idx), error_floor);
        all_values = [all_values, series_values{series_id}(:)']; %#ok<AGROW>
    end

    valid_values = all_values(isfinite(all_values) & all_values > 0);
    if isempty(valid_values)
        valid_values = error_floor;
    end
    y_min = max(error_floor, min(valid_values) * 0.8);
    y_max = max(valid_values) * 1.2;
    if y_max <= y_min
        y_max = y_min * 10;
    end

    video_path = fullfile(pwd, video_name);
    fig_video = figure('Color', 'w', 'Position', [100, 100, 1280, 720], ...
        'Visible', 'off');
    ax_video = axes('Parent', fig_video);
    hold(ax_video, 'on');
    grid(ax_video, 'on');
    box(ax_video, 'on');

    h_lines = gobjects(series_count, 1);
    is_velocity_error_video = strcmp(video_name, 'Experiment3_Velocity_Error_Animation.mp4');
    for series_id = 1:series_count
        current_line_width = line_width;
        current_color = colors(series_id, :);
        if series_id == series_count
            if is_velocity_error_video
                current_line_width = max(line_width * 0.55, 1.2);
            else
                current_line_width = max(line_width + 1.0, 3.0);
            end
        end
        h_lines(series_id) = plot(ax_video, nan, nan, ...
            'Color', current_color, 'LineStyle', line_styles{series_id}, ...
            'LineWidth', current_line_width);
    end

    set(ax_video, 'YScale', 'log', 'FontSize', 15, 'LineWidth', 1.0);
    xlim(ax_video, [1, sample_count]);
    ylim(ax_video, [y_min, y_max]);
    xlabel(ax_video, 'Timestep');
    ylabel(ax_video, y_axis_label);
    % No legend in the error videos to avoid covering dense curves.

    video_writer = VideoWriter(video_path, 'MPEG-4');
    video_writer.FrameRate = video_fps;
    video_writer.Quality = 100;
    open(video_writer);

    x_axis = 1:sample_count;
    frame_index = round(linspace(1, sample_count, total_frames));
    for frame_id = 1:total_frames
        k = frame_index(frame_id);
        for series_id = 1:series_count
            set(h_lines(series_id), 'XData', x_axis(1:k), ...
                'YData', series_values{series_id}(1:k));
        end
        drawnow;
        writeVideo(video_writer, getframe(fig_video));
    end

    close(video_writer);
    close(fig_video);
    disp(['Video saved as: ', video_path]);
end
function output_signal = smooth_pae_input(input_signal, dt, tau)
    if tau <= 0
        output_signal = input_signal;
        return
    end

    alpha = dt / (tau + dt);
    output_signal = zeros(size(input_signal));
    output_signal(:, 1) = input_signal(:, 1);

    for k = 2:size(input_signal, 2)
        output_signal(:, k) = output_signal(:, k - 1) + ...
            alpha * (input_signal(:, k) - output_signal(:, k - 1));
    end
end
function output_signal = limit_vector_norm(input_signal, limit_value)
    if ~isfinite(limit_value) || limit_value <= 0
        output_signal = input_signal;
        return
    end

    output_signal = input_signal;
    for k = 1:size(input_signal, 2)
        vec_norm = norm(input_signal(:, k));
        if vec_norm > limit_value
            output_signal(:, k) = input_signal(:, k) * (limit_value / vec_norm);
        end
    end
end

function velocity = second_order_differentiator(position_signal, dt, tau)
    raw_velocity = zeros(size(position_signal));
    sample_count = size(position_signal, 2);
    if sample_count < 2
        velocity = raw_velocity;
        return
    end

    raw_velocity(:, 1) = (position_signal(:, 2) - position_signal(:, 1)) / dt;
    raw_velocity(:, sample_count) = ...
        (position_signal(:, sample_count) - position_signal(:, sample_count - 1)) / dt;
    if sample_count > 2
        raw_velocity(:, 2:sample_count - 1) = ...
            (position_signal(:, 3:sample_count) - position_signal(:, 1:sample_count - 2)) / (2 * dt);
    end

    alpha = dt / (tau + dt);
    velocity = cascade_second_order_link(raw_velocity, alpha);
    velocity = fliplr(cascade_second_order_link(fliplr(velocity), alpha));
end

function output_signal = cascade_second_order_link(input_signal, alpha)
    first_stage = zeros(size(input_signal));
    output_signal = zeros(size(input_signal));
    first_stage(:, 1) = input_signal(:, 1);
    output_signal(:, 1) = input_signal(:, 1);

    for k = 2:size(input_signal, 2)
        first_stage(:, k) = first_stage(:, k - 1) + ...
            alpha * (input_signal(:, k) - first_stage(:, k - 1));
        output_signal(:, k) = output_signal(:, k - 1) + ...
            alpha * (first_stage(:, k) - output_signal(:, k - 1));
    end
end

function R_WB = rpy_to_R_WB(roll, pitch, yaw)
    cr = cos(roll);
    sr = sin(roll);
    cp = cos(pitch);
    sp = sin(pitch);
    cy = cos(yaw);
    sy = sin(yaw);

    R_WB = [ ...
        cy * cp, cy * sp * sr - sy * cr, cy * sp * cr + sy * sr; ...
        sy * cp, sy * sp * sr + cy * cr, sy * sp * cr - cy * sr; ...
        -sp,     cp * sr,                cp * cr];
end





