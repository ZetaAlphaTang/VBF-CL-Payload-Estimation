function result = VBF_without_VB(observations,cfg)

result=run_filter(observations,cfg,true,false,'CL-pseudolinear');
end
