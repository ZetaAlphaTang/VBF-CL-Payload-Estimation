function result=baseline_filter(obs,cfg,use_pae,kind)
[x,P,Q,F,G]=initial_conditions(obs,cfg);
initial_x=x;initial_P=P;
N=numel(obs.t);X=nan(6,N);
if strcmp(kind,'SHF')
    R0=nominal_covariance(obs,cfg,'CL-pseudolinear',1);
else
    R0=nominal_covariance(obs,cfg,'bearing-only',1);
    if strcmp(kind,'EKF'),R0=R0(1,1)*eye(3);end
end
R=R0;first_indefinite=NaN;failure='';
for k=1:N
    try
        if k==1
            xp=x;Pp=P;
        else
            u=zeros(3,1);if use_pae,u=obs.paeInput(:,k-1);end
            xp=F*x+G*u;Pp=F*P*F'+Q;
        end
        if strcmp(kind,'SHF')
            H=[eye(3),zeros(3)];
            innovation=obs.pUav(:,k)+obs.length(k)*obs.bearing(:,k)-H*xp;
            if k>1
                if cfg.shf.rho==1,weight=1/k;
                else,weight=(1-cfg.shf.rho)/(1-cfg.shf.rho^k);end
                R=(1-weight)*R+weight*(innovation*innovation'-H*Pp*H');
            end
            if isnan(first_indefinite) && min(eig((R+R')/2))<=0
                first_indefinite=k;
            end
        elseif strcmp(kind,'EKF')
            r=xp(1:3)-obs.pUav(:,k);distance=norm(r);
            assert(distance>0 && isfinite(distance),'Undefined predicted range.');
            b=r/distance;H=[(eye(3)-b*b')/distance,zeros(3)];
            innovation=obs.bearing(:,k)-b;
        else
            [z,H]=linearized_measurement(xp,obs,k,'bearing-only');
            innovation=z-H*xp;
        end
        K=Pp*H'/(H*Pp*H'+R);
        x=xp+K*innovation;P=(eye(6)-K*H)*Pp;
        assert(all(isfinite(x)) && all(isfinite(P),'all'),'Nonfinite state/covariance.');
        X(:,k)=x;
    catch error_info
        failure=error_info.message;break;
    end
end
result=struct('X',X,'initialState',initial_x,'initialCovariance',initial_P, ...
    'Q',Q,'Rinitial',R0,'Rfinal',R,'FirstIndefiniteR',first_indefinite,'Failure',failure);
if strcmp(kind,'SHF'),result.ForgettingFactor=cfg.shf.rho;end
end
