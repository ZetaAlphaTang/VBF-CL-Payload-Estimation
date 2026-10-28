function result=EKF(observations,cfg)
result=baseline_filter(observations,cfg,true,'EKF');
end
