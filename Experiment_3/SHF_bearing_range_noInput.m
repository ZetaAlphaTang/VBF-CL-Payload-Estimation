function [est_target_state, P, R] = SHF_bearing_range_noInput(P, state_old, P_t, dt, Q, R, k, rho)
% Sage-Husa pseudo-linear filter without PAE input

w1 = 1 * pi / 180;

F = [1, 0, 0, dt, 0, 0;
    0, 1, 0, 0, sin(w1 * dt) / w1, (cos(w1 * dt) - 1) / w1;
    0, 0, 1, 0, (1 - cos(w1 * dt)) / w1, sin(w1 * dt) / w1;
    0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, cos(w1 * dt), sin(w1 * dt);
    0, 0, 0, 0, -sin(w1 * dt), cos(w1 * dt)];

H = [eye(3), zeros(3)];

% Time update
X_pred = F * state_old;
P_pred = F * P * F' + Q;

% Innovation
innovation = P_t - H * X_pred;

% Sage-Husa adaptive measurement covariance
beta = (1 - rho) / (1 - rho^k);
R_new = innovation * innovation' - H * P_pred * H';
R = (1 - beta) * R + beta * R_new;
R = regularize_cov(R, 1e-8);

% Measurement update
S = H * P_pred * H' + R;
K = P_pred * H' * pinv(S);

est_target_state = X_pred + K * innovation;

I6 = eye(6);
P = (eye(6) - K * H) * P_pred;
end