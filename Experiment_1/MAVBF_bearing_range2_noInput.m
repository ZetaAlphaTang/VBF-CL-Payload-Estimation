function [Xkk, Pkk, jkk, Jkk, e, S] = MAVBF_bearing_range2_noInput(Pkk, P_t, l, Xkk, dt, Q, jkk, Jkk, bearing)
% Variational Bayesian pseudo-linear filter without PAE input

n = 6;
m = 3;

w1 = 1 * pi / 180;

F = [1, 0, 0, dt, 0, 0;
    0, 1, 0, 0, sin(w1 * dt) / w1, (cos(w1 * dt) - 1) / w1;
    0, 0, 1, 0, (1 - cos(w1 * dt)) / w1, sin(w1 * dt) / w1;
    0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, cos(w1 * dt), sin(w1 * dt);
    0, 0, 0, 0, -sin(w1 * dt), cos(w1 * dt)];

H = [eye(3), zeros(3)];

% Nominal time update without acceleration input
Xkk_1 = F * Xkk;
Pkk_1 = F * Pkk * F' + Q;

% VB initialization
Xkk = Xkk_1;
Pkk = Pkk_1;

rho = 1 - exp(-4);
alpha = 3.6;

dkk_1 = n + alpha + 1;
Dkk_1 = alpha * Pkk_1;

jkk_1 = rho * (jkk - m - 1) + m + 1;
Jkk_1 = rho * Jkk;

for Ll = 1 : 1

    % Update predicted covariance P(k+1|k)
    Mkk = (Xkk - Xkk_1) * (Xkk - Xkk_1)' + Pkk;

    dkk = dkk_1 + 1;
    Dkk = Dkk_1 + Mkk;

    % Update measurement covariance R(k+1)
    innovation_tmp = P_t - H * Xkk;
    Okk = innovation_tmp * innovation_tmp' + H * Pkk * H';

    jkk = jkk_1 + 1;
    Jkk = Jkk_1 + Okk;

    % Effective covariance matrices
    Pkk_1_eff = Dkk / (dkk - n - 1);
    R_eff = Jkk / (jkk - m - 1);

    Pkk_1_eff = regularize_cov(Pkk_1_eff, 1e-10);
    R_eff = regularize_cov(R_eff, 1e-8);

    % State update
    S = H * Pkk_1_eff * H' + R_eff;
    K = Pkk_1_eff * H' * pinv(S);

    e = P_t - H * Xkk_1;
    Xkk = Xkk_1 + K * e;

    I6 = eye(6);
    Pkk = (I6 - K * H) * Pkk_1_eff * (I6 - K * H)' + K * R_eff * K';
    Pkk = regularize_cov(Pkk, 1e-10);
end
end