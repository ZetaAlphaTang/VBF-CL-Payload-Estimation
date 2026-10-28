function [X_PCMPC_EKF, parameters] = PCMPC_EKF(bearing, P_observation, P_UAV, attitude_UAV, T_pred, T, l, m_quad, estimation_idx)
cfg = struct('bearing_sigma',1.5 * pi/180,'acceleration_sigma',0.5);
P0 =  1 * diag([(1.0*pi/180)^2*ones(1,3),(1.0*pi/180)^2*ones(1,3)]);
G = [0.5*T^2*eye(3);T*eye(3)];
parameters = cfg;
parameters.P0 = P0;
parameters.Q = G*((cfg.acceleration_sigma/l)^2*eye(3))*G';
position_variance = (l*cfg.bearing_sigma)^2;
parameters.R = blkdiag(position_variance*eye(3),(2*position_variance/T^2)*eye(3));
parameters.derivative_tau = 0;
baseline_nav_velocity = UAV_velocity(P_UAV,T);
baseline_total_samples = size(P_UAV,2);
baseline_first_sample = estimation_idx(1);
baseline_R_WB = zeros(3,3,baseline_total_samples);
for i = 1:baseline_total_samples
    baseline_R_WB(:,:,i) = Attitude_rotation(attitude_UAV(1,i),attitude_UAV(2,i),attitude_UAV(3,i));
end
X_PCMPC_internal = [bearing(:,baseline_first_sample); zeros(3,1)];
P_PCMPC = P0;
X_PCMPC_EKF = nan(6, baseline_total_samples);
pcm_position_body_previous = l * baseline_R_WB(:,:,baseline_first_sample)' * bearing(:,baseline_first_sample);
for i = estimation_idx
    step_dt = T * (i ~= baseline_first_sample);
    input_idx = max(baseline_first_sample, i - 1);
    pcm_position_body = l * baseline_R_WB(:,:,i)' * bearing(:,i);
    if i == baseline_first_sample
        pcm_velocity_body_measurement = zeros(3,1);
        pcm_omega_body = zeros(3,1);
    else
        pcm_velocity_body_measurement = (pcm_position_body-pcm_position_body_previous)/T;
        pcm_omega_body = UAV_angular_velocity(baseline_R_WB(:,:,i-1),baseline_R_WB(:,:,i),T);
    end
    pcm_position_body_previous = pcm_position_body;
    pcm_thrust_world = T_pred(input_idx)*baseline_R_WB(:,3,input_idx);
    [X_PCMPC_internal, P_PCMPC] = EKF_update( ...
        X_PCMPC_internal, P_PCMPC, pcm_position_body, pcm_velocity_body_measurement, ...
        pcm_thrust_world, baseline_R_WB(:,:,i), pcm_omega_body, step_dt, T, l, m_quad, cfg);
    assert(all(isfinite(X_PCMPC_internal)) && all(isfinite(P_PCMPC),'all'), ...
        'PCMPC_EKF:NonfiniteResult','PCMPC has a nonfinite result at sample %d.',i);
    X_PCMPC_EKF(:,i) = [P_observation(:,i)+l*X_PCMPC_internal(1:3); ...
        baseline_nav_velocity(:,i)+l*X_PCMPC_internal(4:6)];
end
end

function [state_new, covariance_new] = EKF_update(state, covariance, ...
    position_body, velocity_body, thrust_W, R_WB, omega_body, dt, sample_dt, ...
    cable_length, uav_mass, cfg)
    [predicted_state,A] = Euler_prediction(state,thrust_W,dt,cable_length,uav_mass);
    G = [0.5*dt^2*eye(3);dt*eye(3)];
    Q = G*((cfg.acceleration_sigma/cable_length)^2*eye(3))*G';
    predicted_covariance = A*covariance*A'+Q;
    H = [cable_length*R_WB',zeros(3); ...
        -cable_length*Skew_symmetric_matrix(omega_body)*R_WB',cable_length*R_WB'];
    position_variance = (cable_length*cfg.bearing_sigma)^2;
    R_measurement = blkdiag(position_variance*eye(3), ...
        (2*position_variance/sample_dt^2)*eye(3));
    measurement = [position_body;velocity_body];
    if dt == 0
        H = H(1:3,:);
        measurement = position_body;
        R_measurement = R_measurement(1:3,1:3);
    end
    innovation_covariance = H*predicted_covariance*H'+R_measurement;
    gain = predicted_covariance*H'/innovation_covariance;
    state_new = predicted_state+gain*(measurement-H*predicted_state);
    covariance_new = (eye(6)-gain*H)*predicted_covariance;
end

function [predicted_state,A] = Euler_prediction(state,thrust_W,dt,cable_length,uav_mass)
    xi = state(1:3);
    xi_dot = state(4:6);
    xi_ddot = cross(xi,cross(xi,thrust_W))/(uav_mass*cable_length) ...
        - dot(xi_dot,xi_dot)*xi;
    predicted_state = state+dt*[xi_dot;xi_ddot];
    A_xi = (dot(xi,thrust_W)*eye(3)+xi*thrust_W'-2*thrust_W*xi') ...
        /(uav_mass*cable_length)-dot(xi_dot,xi_dot)*eye(3);
    A_rate = -2*xi*xi_dot';
    A = eye(6)+dt*[zeros(3),eye(3);A_xi,A_rate];
end

function omega = UAV_angular_velocity(R_previous,R_current,dt)
    increment = R_previous'*R_current;
    angle = acos(max(-1,min(1,(trace(increment)-1)/2)));
    vector = [increment(3,2)-increment(2,3); ...
        increment(1,3)-increment(3,1);increment(2,1)-increment(1,2)];
    if angle < 1e-7
        omega = vector/(2*dt);
    else
        assert(angle < pi-1e-4,'Attitude increment is too large for gyro reconstruction.');
        omega = angle*vector/(2*sin(angle)*dt);
    end
end

function R_WB = Attitude_rotation(roll,pitch,yaw)
    cr=cos(roll); sr=sin(roll);
    cp=cos(pitch); sp=sin(pitch);
    cy=cos(yaw); sy=sin(yaw);
    R_WB = [cy*cp,cy*sp*sr-sy*cr,cy*sp*cr+sy*sr; ...
        sy*cp,sy*sp*sr+cy*cr,sy*sp*cr-cy*sr;-sp,cp*sr,cp*cr];
end

function velocity = UAV_velocity(position,dt)
    velocity = zeros(size(position));
    velocity(:,2:end) = diff(position,1,2)/dt;
end

function S = Skew_symmetric_matrix(vector)
    S = [0,-vector(3),vector(2);vector(3),0,-vector(1);-vector(2),vector(1),0];
end
