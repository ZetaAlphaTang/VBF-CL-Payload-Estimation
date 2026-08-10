function [est_target_state, P] = EKF_bearing_range_noInput(P, state_old, bearing, P_UAV, dt, Q, R)
% EKF bearing-only estimator without PAE input
% Measurement equation:
% z(k) = g_m(k)
% h(x) = (p_l(k) - p_u(k)) / ||p_l(k) - p_u(k)||

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

% Time update
X_pred = F * state_old;
P_pred = F * P * F' + Q;

% Normalize measured bearing
z = bearing;

% Nonlinear bearing prediction
rel = X_pred(1:3) - P_UAV;
r = norm(rel);

g_pred = rel / r;

% Measurement function
h = g_pred;

% Jacobian of bearing with respect to payload position
H_g = (eye(3) - g_pred * g_pred') / r;

% Full measurement Jacobian with respect to [p_l; v_l]
H = [H_g, zeros(3)];

% Innovation
innovation = z - h;

% Measurement update
S = H * P_pred * H' + R;
K = P_pred * H' * pinv(S);

est_target_state = X_pred + K * innovation;

% Joseph-form covariance update
I6 = eye(6);
P = (eye(6) - K * H) * P_pred;
end