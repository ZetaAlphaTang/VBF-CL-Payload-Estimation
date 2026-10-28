function result = VBF_without_PAE(observations,cfg)

result=run_filter(observations,cfg,false,true,'CL-pseudolinear');
end
