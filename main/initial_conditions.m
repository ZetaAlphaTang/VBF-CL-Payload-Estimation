function [x,P,Q,F,G] = initial_conditions(obs,cfg)
[F,G]=CA_model(cfg.dt);
Q=cfg.vb.qScale*G*(cfg.modelAccelerationSigma^2*eye(3))*G';
b0=obs.bearing(:,1);tangent=eye(3)-b0*b0';
x=[obs.pUav(:,1)+cfg.initialRange*b0;zeros(3,1)];
P=blkdiag(cfg.initialRangeSigma^2*(b0*b0')+ ...
    (cfg.nominal.initialPositionSigma^2+(cfg.initialRange*cfg.nominal.initialBearingSigma)^2)*tangent, ...
    cfg.initialVelocitySigma^2*eye(3));
end
