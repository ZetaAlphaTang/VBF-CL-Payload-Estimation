function result = VBF_without_CL(observations,cfg)

result=run_filter(observations,cfg,true,true,'bearing-only');
end
