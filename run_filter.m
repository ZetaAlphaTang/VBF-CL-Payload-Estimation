function result = run_filter(obs,cfg,use_pae,use_vb,model)

[x,P,Q,F,G] = initial_conditions(obs,cfg);
initial_x=x; initial_P=P;
if strcmp(model,'bearing-only'), m=2; else, m=3; end
R0=nominal_covariance(obs,cfg,model,1);
j=m+1+cfg.vb.equivalentSamples;
J=cfg.vb.equivalentSamples*R0;
N=numel(obs.t); X=zeros(6,N);
for k=1:N
    if k==1
        xp=x; Pp=P;
    else
        u=zeros(3,1);
        if use_pae, u=obs.paeInput(:,k-1); end
        xp=F*x+G*u;
        Pp=regularize_cov(F*P*F'+Q,cfg.epsilon);
    end
    [z,H]=linearized_measurement(xp,obs,k,model);
    if use_vb
        [x,P,j,J,R]=vb_update(xp,Pp,z,H,j,J,cfg);
    else

        R=nominal_covariance(obs,cfg,model,k);
        [x,P]=gaussian_update(xp,Pp,z,H,R,cfg.epsilon);
    end
    X(:,k)=x;
    if any(~isfinite(x)) || any(~isfinite(P),'all')
        X(:,k:end)=NaN;
        break;
    end
end
result=struct('X',X, ...
    'initialState',initial_x,'initialCovariance',initial_P, ...
    'Rinitial',R0,'Q',Q);
end
