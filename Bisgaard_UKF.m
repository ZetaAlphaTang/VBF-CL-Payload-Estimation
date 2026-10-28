function [X_Bisgaard_UKF, parameters] = Bisgaard_UKF(bearing, P_observation, P_UAV, IMU_UAV, T, l, estimation_idx)
cfg = struct('bearing_sigma',1.5*pi/180,'acceleration_sigma',0.5, ...
    'alpha',0.5,'beta',0.0,'kappa',1.0,'gravity',9.81);
P0 = 0.01 * diag([(1.0*pi/180)^2*ones(1,2),(1.0*pi/180)^2*ones(1,2)]);
baseline_total_samples = size(P_UAV,2);
baseline_first_sample = estimation_idx(1);
if isscalar(l), l = repmat(l,1,baseline_total_samples); end
l = reshape(l,1,[]);
assert(numel(l)==baseline_total_samples && all(isfinite(l)) && all(l>0), ...
    'Bisgaard_UKF:InvalidLength','Cable-length inputs must be finite and positive.');
G = [0.5*T^2*eye(2);T*eye(2)];
parameters = cfg;
parameters.P0 = P0;
parameters.Q = 10 * G*((cfg.acceleration_sigma/l(baseline_first_sample))^2*eye(2))*G';
parameters.R = cfg.bearing_sigma^2*eye(3);
parameters.derivative_tau = 0;
parameters.CableLength_m = l;
parameters.Q_used = nan(4,4,baseline_total_samples);
baseline_nav_velocity = UAV_velocity(P_UAV,T);
initial_down_bearing = diag([1,1,-1]) * bearing(:,baseline_first_sample);
X_Bisgaard_internal = [atan2(initial_down_bearing(1),initial_down_bearing(3)); ...
    asin(initial_down_bearing(2));0;0];
P_Bisgaard = P0;
X_Bisgaard_UKF = nan(6,baseline_total_samples);
for i = estimation_idx
    step_dt = T * (i ~= baseline_first_sample);
    input_idx = max(baseline_first_sample,i-1);
    [X_Bisgaard_internal,P_Bisgaard] = UKF_update( ...
        X_Bisgaard_internal,P_Bisgaard,bearing(:,i), ...
        IMU_UAV(:,input_idx),step_dt,l(input_idx),cfg);
    [relative_position,relative_velocity] = Payload_state(X_Bisgaard_internal,l(i));
    G_step = [0.5*step_dt^2*eye(2);step_dt*eye(2)];
    parameters.Q_used(:,:,i) = 10*G_step*((cfg.acceleration_sigma/l(input_idx))^2*eye(2))*G_step';
    X_Bisgaard_UKF(:,i) = [P_observation(:,i)+relative_position; ...
        baseline_nav_velocity(:,i)+relative_velocity];
    assert(all(isfinite(X_Bisgaard_UKF(:,i))), ...
        'Bisgaard_UKF:NonfiniteState','Nonfinite state at sample %d.',i);
end
end

function [state_new,covariance_new] = UKF_update( ...
    state,covariance,measured_bearing_W,uav_acceleration_W,dt,cable_length,cfg)
n = 4;
lambda = cfg.alpha^2*(n+cfg.kappa)-n;
mean_weights = [lambda/(n+lambda),repmat(1/(2*(n+lambda)),1,2*n)];
covariance_weights = mean_weights;
covariance_weights(1) = covariance_weights(1)+1-cfg.alpha^2+cfg.beta;
sigma = Sigma_points(state,covariance,n+lambda);
for sigma_id = 1:size(sigma,2)
    sigma(:,sigma_id) = sigma(:,sigma_id)+dt*Pendulum_model( ...
        sigma(:,sigma_id),uav_acceleration_W,cable_length,cfg.gravity);
end
predicted_state = sigma*mean_weights';
G = [0.5*dt^2*eye(2);dt*eye(2)];
predicted_covariance = 10*G*((cfg.acceleration_sigma/cable_length)^2*eye(2))*G';
for sigma_id = 1:size(sigma,2)
    delta = sigma(:,sigma_id)-predicted_state;
    predicted_covariance = predicted_covariance+covariance_weights(sigma_id)*(delta*delta');
end
sigma = Sigma_points(predicted_state,predicted_covariance,n+lambda);
predicted_measurements = zeros(3,size(sigma,2));
for sigma_id = 1:size(sigma,2)
    predicted_measurements(:,sigma_id) = Cable_direction(sigma(:,sigma_id));
end
measurement_mean = predicted_measurements*mean_weights';
S = cfg.bearing_sigma^2*eye(3);
C = zeros(4,3);
for sigma_id = 1:size(sigma,2)
    dx = sigma(:,sigma_id)-predicted_state;
    dz = predicted_measurements(:,sigma_id)-measurement_mean;
    S = S+covariance_weights(sigma_id)*(dz*dz');
    C = C+covariance_weights(sigma_id)*(dx*dz');
end
gain = C/S;
state_new = predicted_state+gain*(measured_bearing_W-measurement_mean);
covariance_new = predicted_covariance-gain*S*gain';
end

function derivative = Pendulum_model(state,acceleration_W,cable_length,gravity)
acceleration_down = diag([1,1,-1])*acceleration_W;
theta = state(1);
phi = state(2);
ax = acceleration_down(1);
ay = acceleration_down(2);
az = acceleration_down(3);
theta_ddot = (-cos(theta)*cos(phi)*ax+ ...
    sin(theta)*cos(phi)*(az-gravity))/cable_length;
phi_ddot = (sin(theta)*sin(phi)*ax-cos(phi)*ay+ ...
    cos(theta)*sin(phi)*(az-gravity))/cable_length;
derivative = [state(3:4);theta_ddot;phi_ddot];
end

function bearing_W = Cable_direction(state)
theta = state(1);
phi = state(2);
bearing_W = [sin(theta)*cos(phi);sin(phi);-cos(theta)*cos(phi)];
end

function [relative_position,relative_velocity] = Payload_state(state,cable_length)
theta = state(1);
phi = state(2);
J = [cos(theta)*cos(phi),-sin(theta)*sin(phi); ...
    0,cos(phi);sin(theta)*cos(phi),cos(theta)*sin(phi)];
relative_position = cable_length*Cable_direction(state);
relative_velocity = cable_length*J*state(3:4);
end

function points = Sigma_points(mean_state,covariance,spread)
root = chol(spread*covariance,'lower');
points = [mean_state,mean_state+root,mean_state-root];
end

function velocity = UAV_velocity(position,dt)
velocity = [zeros(3,1),diff(position,1,2)/dt];
end
