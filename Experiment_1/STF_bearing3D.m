function [est_target_state, P, V_k] = STF_bearing3D(P, bearing, state_old, P_t, dt, Q, R, k, rou, V_k, u,l)
    % bearing-only
    % pseudo linear kalman filter estimator
    Xk_1 = state_old;
    Pk_1 = P;
    % state transition matrix
    w1 = 1 * pi / 180; % CT angular-rate parameter
    F = [1, 0, 0, dt, 0, 0;
        0, 1, 0, 0, sin(w1 * dt) / w1, (cos(w1 * dt) - 1)/ w1;
        0, 0, 1, 0, (1 - cos(w1 * dt))/ w1, sin(w1 * dt) / w1;
        0, 0, 0, 1, 0, 0;
        0, 0, 0, 0, cos(w1 * dt), sin(w1 * dt);
        0, 0, 0, 0, -sin(w1 * dt), cos(w1 * dt)];
    % pseudo-linear measurement matrix
    G = [0.5 * dt *dt, 0, 0; 0, 0.5 * dt *dt, 0; 0, 0, 0.5 * dt *dt; dt, 0, 0; 0, dt, 0; 0, 0, dt];
    % pseudo-linear measurement matrix
    H = [eye(3), zeros(3); zeros(3), zeros(3)];
    % pseudolinear measurements
    mear_state = [P_t; u];
    % noise covariance
    G_H = [zeros(3); eye(3)];
    E = [l * eye(3), bearing, zeros(3); zeros(3), zeros(3,1), eye(3)];
    R = E * R * E';
    % Time Update
    Xkk_1 = F * Xk_1  + G * u;
    Zkk_1 = H * Xkk_1;
    gamma = mear_state - Zkk_1 - G_H * u;
    if k == 1
        V_k = gamma * gamma';
    else
        V_k = (rou * V_k + gamma * gamma') / (1 + rou);
    end
    N_k = V_k - H * Q * H' - R;
    M_k = H * F * Pk_1 * F' * H';
    lam_0 = trace(N_k) / trace(M_k);
    
    if lam_0 >= 1
        lam_k = lam_0;
    else
        lam_k = 1;
    end

    Pkk_1 = lam_k * F * P * F' + Q;
    % Kalman filter steps
    K = Pkk_1 * H' * pinv(H * Pkk_1 * H' + R);
    P = (eye(6) - K * H) * Pkk_1;
    est_target_state = Xkk_1 + K * gamma;

end

