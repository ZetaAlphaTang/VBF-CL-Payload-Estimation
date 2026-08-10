function [est_target_state, P, R] = SHF_bearing_range(P, state_old, P_t, dt, Q, R, u, k, rho)
% Sage-Husa pseudo-linear filter with PAE input

w1 = 1 * pi / 180;

F = [1, 0, 0, dt, 0, 0;
    0, 1, 0, 0, sin(w1 * dt) / w1, (cos(w1 * dt) - 1) / w1;
    0, 0, 1, 0, (1 - cos(w1 * dt)) / w1, sin(w1 * dt) / w1;
    0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, cos(w1 * dt), sin(w1 * dt);
    0, 0, 0, 0, -sin(w1 * dt), cos(w1 * dt)];

G = [0.5 * dt * dt, 0, 0;
    0, 0.5 * dt * dt, 0;
    0, 0, 0.5 * dt * dt;
    dt, 0, 0;
    0, dt, 0;
    0, 0, dt];

H = [eye(3), zeros(3)];

% Time update
X_pred = F * state_old + G * u;
P_pred = F * P * F' + Q;

% Innovation
innovation = P_t - H * X_pred;

% Sage-Husa adaptive measurement covariance
beta = (1 - rho) / (1 - rho^k);
R_new = innovation * innovation' - H * P_pred * H';
R = (1 - beta) * R + beta * R_new;
R = bound_cov(R, 1e-8, 0.005);

% Measurement update
S = H * P_pred * H' + R;
K = P_pred * H' * pinv(S);

est_target_state = X_pred + K * innovation;

I6 = eye(6);
P = (I6 - K * H) * P_pred * (I6 - K * H)' + K * R * K';
P = regularize_cov(P, 1e-12);
end

function C = bound_cov(C, lower_value, upper_value)
    C = (C + C') / 2;
    [V, D] = eig(C);
    d = max(real(diag(D)), lower_value);
    if isfinite(upper_value)
        d = min(d, upper_value);
    end
    C = V * diag(d) * V';
    C = (C + C') / 2;
end


