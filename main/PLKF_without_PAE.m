function result=PLKF_without_PAE(observations,cfg)
result=run_filter(observations,cfg,false,false,'CL-pseudolinear');
end
