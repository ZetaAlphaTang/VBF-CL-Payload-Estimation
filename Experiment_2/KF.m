function [est_target_state, P] = KF(P, state_old, P_t, dt,Q,R,u)
% bearing-range
% pseudo linear kalman filter estimator
w1 = 1 * pi / 180; % CT angular-rate parameter
% state transition matrix
F = [1, 0, 0, dt, 0, 0;
    0, 1, 0, 0, sin(w1 * dt) / w1, (cos(w1 * dt) - 1)/ w1;
    0, 0, 1, 0, (1 - cos(w1 * dt))/ w1, sin(w1 * dt) / w1;
    0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, cos(w1 * dt), sin(w1 * dt);
    0, 0, 0, 0, -sin(w1 * dt), cos(w1 * dt)];
G = [0.5 * dt *dt, 0, 0; 0, 0.5 * dt *dt, 0; 0, 0, 0.5 * dt *dt; dt, 0, 0; 0, dt, 0; 0, 0, dt];
% pseudo-linear measurement matrix
H = [eye(3), zeros(3)];

% pseudolinear measurements
mear_state = P_t;
% noise covariance
% Kalman filter steps
state_old = F * state_old + G * u;
P = F * P * F' + Q;
K = P * H' * pinv(H * P * H' + R);
est_target_state = state_old + K * (mear_state - H * state_old);
P = (eye(6) - K * H) * P;
end

