function result = VBF_CL(observations,cfg)

result=run_filter(observations,cfg,true,true,'CL-pseudolinear');
end
