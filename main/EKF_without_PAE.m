function result=EKF_without_PAE(observations,cfg)
result=baseline_filter(observations,cfg,false,'EKF');
end
