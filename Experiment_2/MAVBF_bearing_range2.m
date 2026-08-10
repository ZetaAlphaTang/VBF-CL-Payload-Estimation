function [Xkk, Pkk, jkk, Jkk, e, S] = MAVBF_bearing_range2(Pkk, P_t, l, Xkk, dt, Q, jkk, Jkk, u, bearing)
% bearing-only
% pseudo linear kalman filter estimator
n = 6;
m = 3;
% state transition matrix
w1 = 1 * pi / 180; % CT angular-rate parameter
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
G_H = [zeros(3)];
E = [l * eye(3), bearing, zeros(3)];

% Kalman filter steps
Xkk_1 = F * Xkk + G * u;
Pkk_1 = F * Pkk * F' + Q;

%% VB initialization
Xkk = Xkk_1;
Pkk = Pkk_1;
rou = 1 - exp(-0.1);% Forgetting factor
tao1 = 3.6;
ikk_1 = (n + 1 + tao1);
Ikk_1 = tao1 * Pkk_1;
jkk_1 = rou * (jkk - m - 1) + m + 1;
Jkk_1 = rou * Jkk;
for Ll= 1 : 3
    % Update predicted covariance P(k+1|k)
    Mkk = (Xkk_1 - Xkk) * (Xkk_1 - Xkk)' + Pkk;% Predicted covariance residual
    ikk = ikk_1 + 1;
    Ikk = Ikk_1 + Mkk;
    % Update measurement covariance R(k+1)
    Okk = (mear_state - H * Xkk) * (mear_state - H * Xkk)' + H * Pkk *H';%(E' * pinv(E * E') * mear_state - E' * pinv(E * E') * H * Xkk) * (E' * pinv(E * E') * mear_state - E' * pinv(E * E') * H * Xkk)' + (E' * pinv(E * E') * H)*Pkk*(E' * pinv(E * E') * H)';%R
    jkk = jkk_1 + 1;
    Jkk = Jkk_1 + Okk;
    % Compute effective P(k+1|k) and R(k+1)
    Pkk_1 = regularize_cov((Ikk) / (ikk - n - 1), 1e-12);
    R = regularize_cov((Jkk) / (jkk - m - 1), 1e-12);
    % Update state and covariance
    S = H * Pkk_1 * H' + R;  % Innovation covariance
    %S = S + 1e-2 * eye(size(S));
    K = Pkk_1 * H' * pinv(S);
    e = mear_state - H * Xkk_1;
    Xkk = Xkk_1 + K * e;
    
    I6 = eye(6);
    Pkk = (I6 - K * H) * Pkk_1 * (I6 - K * H)' + K * R * K';
    Pkk = regularize_cov(Pkk, 1e-12);
end
end



