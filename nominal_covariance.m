function R=nominal_covariance(obs,cfg,model,k)

sp2=cfg.nominal.uavPositionVariance;
sc2=cfg.nominal.cableVariance;
sb2=cfg.nominal.angularVariance;
if strcmp(model,'bearing-only')

    R=(sb2+sp2/cfg.initialRange^2)*eye(2);
elseif strcmp(model,'CL-pseudolinear')
    b=obs.bearing(:,k); l=cfg.cableLength;
    R=sp2*eye(3)+sc2*(b*b')+l^2*sb2*(eye(3)-b*b');
else
    error('Unknown measurement model.');
end

R=regularize_cov(cfg.vb.rScale*R,cfg.epsilon);
end
