function data = generate_simulation(cfg)
t=cfg.t; dt=cfg.dt; N=numel(t); l=cfg.cableLength; g=cfg.g;
[pU,vU,aU,q0,qd0]=figure_eight_motion(t);
q=zeros(3,N); qd=q; qdd=q;
q(:,1)=q0/norm(q0); qd(:,1)=qd0-q(:,1)*(q(:,1)'*qd0);

for k=1:N
    qdd(:,k)=pendulum_acceleration(q(:,k),qd(:,k),aU(:,k),l,g);
    if k<N
        y=[q(:,k);qd(:,k)]; h=dt/2;
        for sub=1:2
            s=(sub-1)/2;
            a0=(1-s)*aU(:,k)+s*aU(:,k+1);
            am=(1-s-0.25)*aU(:,k)+(s+0.25)*aU(:,k+1);
            a1=(1-s-0.5)*aU(:,k)+(s+0.5)*aU(:,k+1);
            k1=pendulum_rhs(y,a0,l,g); k2=pendulum_rhs(y+h*k1/2,am,l,g);
            k3=pendulum_rhs(y+h*k2/2,am,l,g); k4=pendulum_rhs(y+h*k3,a1,l,g);
            y=y+h*(k1+2*k2+2*k3+k4)/6;
            y(1:3)=y(1:3)/norm(y(1:3));
            y(4:6)=y(4:6)-y(1:3)*(y(1:3)'*y(4:6));
        end
        q(:,k+1)=y(1:3); qd(:,k+1)=y(4:6);
    end
end
pL=pU+l*q; vL=vU+l*qd; aL=aU+l*qdd;
fCable=cfg.mPayload*(g-aL);
thrustVector=cfg.mUav*(aU-g)-fCable;
thrust=vecnorm(thrustVector); bz=thrustVector./thrust;
assert(all(sum(fCable.*q,1)>0),'Truth has slack/negative-tension cable.');

pos_normal=randn(3,N); angular_normal=randn(3,N); length_normal=randn(1,N);
imu_normal=randn(3,N); thrust_normal=randn(1,N); attitude_normal=randn(3,N);
outlier_uniform=rand(1,N);
angular_state=zeros(3,N); angular_state(:,1)=angular_normal(:,1);
for k=2:N
    angular_state(:,k)=cfg.sensor.bearingCorrelation*angular_state(:,k-1)+ ...
        sqrt(1-cfg.sensor.bearingCorrelation^2)*angular_normal(:,k);
end
envelope=0.5+2.5*sin(2*pi*0.12*t+0.4).^2;
visual_envelope=cfg.sensor.visualEnvelopeMinimum+ ...
    cfg.sensor.visualEnvelopeAmplitude*sin(2*pi*0.12*t+0.4).^2;
visual_scale=visual_envelope;
visual_scale(outlier_uniform<cfg.sensor.outlierProbability)= ...
    visual_scale(outlier_uniform<cfg.sensor.outlierProbability)*cfg.sensor.outlierScale;
pMeasured=pU+cfg.sensor.uavPositionSigma*pos_normal.*envelope;
bMeasured=zeros(3,N); bzMeasured=zeros(3,N); specificForceMeasured=zeros(3,N);
attitudeMeasured=zeros(3,N); rotationMeasured=zeros(3,3,N);
for k=1:N
    deltaAttitude=cfg.sensor.attitudeSigma*attitude_normal(:,k);
    visualRotationError=cfg.sensor.bearingBias+ ...
        cfg.sensor.bearingSigma*visual_scale(k)* ...
        cfg.sensor.bearingNoiseShape*angular_state(:,k);
    cameraBearing=rotate_vector(q(:,k),visualRotationError);
    bMeasured(:,k)=rotate_vector(cameraBearing,deltaAttitude);
    bzMeasured(:,k)=rotate_vector(bz(:,k),deltaAttitude);

    specificForceMeasured(:,k)=rotate_vector(aU(:,k)-g+ ...
        cfg.sensor.imuSigma*imu_normal(:,k),deltaAttitude);
    Rtrue=thrust_attitude(bz(:,k));
    Rmeasured=[rotate_vector(Rtrue(:,1),deltaAttitude), ...
        rotate_vector(Rtrue(:,2),deltaAttitude),bzMeasured(:,k)];
    rotationMeasured(:,:,k)=Rmeasured;
    attitudeMeasured(:,k)=[atan2(Rmeasured(3,2),Rmeasured(3,3)); ...
        asin(max(-1,min(1,-Rmeasured(3,1)))); ...
        atan2(Rmeasured(2,1),Rmeasured(1,1))];
end
lMeasured=l+cfg.sensor.cableSigma*length_normal.*envelope;
assert(all(lMeasured>0),'Nonpositive length measurement; do not clip it silently.');

thrustMeasured=thrust.*(1+cfg.sensor.thrustRelativeSigma*thrust_normal);

pae=PAE(specificForceMeasured,thrustMeasured,bzMeasured,cfg.mUav,cfg.mPayload,g);
paeNoiseless=PAE(aU-g,thrust,bz,cfg.mUav,cfg.mPayload,g);
E=zeros(3,2,N); E(:,:,1)=tangent_basis(bMeasured(:,1));
for k=2:N

    E(:,:,k)=transport_basis(E(:,:,k-1),bMeasured(:,k-1),bMeasured(:,k));
end
obs=struct('t',t,'pUav',pMeasured,'bearing',bMeasured,'length',lMeasured, ...
    'paeInput',pae,'tangent',E,'accelerationUav',specificForceMeasured+g, ...
    'specificForceWorld',specificForceMeasured,'thrust',thrustMeasured, ...
    'attitude',attitudeMeasured,'rotationWorldBody',rotationMeasured);

data=struct('t',t,'pPayloadTrue',pL,'vPayloadTrue',vL,'aPayloadTrue',aL, ...
    'pUavTrue',pU,'vUavTrue',vU,'aUavTrue',aU,'observations',obs, ...
    'visualScale',visual_scale,'noiseEnvelope',envelope, ...
    'noiselessPaeResidual',max(vecnorm(paeNoiseless-aL)), ...
    'cableConstraintResidual',max(abs(vecnorm(pL-pU)-l)), ...
    'specificForceMeasured',specificForceMeasured,'thrustMeasured',thrustMeasured, ...
    'bodyZTrue',bz,'bodyZMeasured',bzMeasured,'cableDirectionTrue',q);
assert(data.noiselessPaeResidual<1e-9,'PAE sign/gravity/force convention failed.');
assert(data.cableConstraintResidual<1e-10,'Cable constraint failed.');
end

function R=thrust_attitude(bz)

roll=-asin(max(-1,min(1,bz(2)))); pitch=atan2(bz(1),bz(3));
cr=cos(roll); sr=sin(roll); cp=cos(pitch); sp=sin(pitch);
R=[cp,sp*sr,sp*cr;0,cr,-sr;-sp,cp*sr,cp*cr];
end

function dy=pendulum_rhs(y,aU,l,g)
dy=[y(4:6);pendulum_acceleration(y(1:3),y(4:6),aU,l,g)];
end

function qdd=pendulum_acceleration(q,qd,aU,l,g)
qdd=(eye(3)-q*q')*(g-aU)/l-(qd'*qd)*q;
end

function [p,v,a,q0,qd0]=figure_eight_motion(t)
w=2*pi/10;
p=[.48*sin(w*t);.30*sin(2*w*t);1.22+.02*sin(.5*w*t)];
v=[.48*w*cos(w*t);.60*w*cos(2*w*t);.01*w*cos(.5*w*t)];
a=[-.48*w^2*sin(w*t);-1.20*w^2*sin(2*w*t);-.005*w^2*sin(.5*w*t)];
q0=[.14;-.12;-1]; qd0=[.38;.28;0];
end

function y=rotate_vector(x,d)
a=norm(d);
if a<1e-14, y=x; return; end
u=d/a; y=x*cos(a)+cross(u,x)*sin(a)+u*(u'*x)*(1-cos(a));

end

function E=tangent_basis(b)
[~,j]=min(abs(b)); e=zeros(3,1); e(j)=1;
e=e-b*(b'*e); e=e/norm(e); E=[e,cross(b,e)];
end

function E=transport_basis(E,bOld,bNew)
v=cross(bOld,bNew); a=atan2(norm(v),dot(bOld,bNew));
if norm(v)>1e-12
    E=[rotate_vector(E(:,1),a*v/norm(v)),rotate_vector(E(:,2),a*v/norm(v))];
elseif dot(bOld,bNew)<0
    E=tangent_basis(bNew);
end

e=E(:,1)-bNew*(bNew'*E(:,1)); e=e/norm(e); E=[e,cross(bNew,e)];
end
