function result = VBF_CL_uncertain_length(observations,cfg,length_bias_m)

assert(isscalar(length_bias_m) && isfinite(length_bias_m));
observations.length=observations.length+length_bias_m;
assert(all(observations.length>0),'Length perturbation is nonpositive; no clipping.');
result=VBF_CL(observations,cfg);
end
