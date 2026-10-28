function [x,P,j,J,R]=vb_update(xp,Pp,z,H,jPrevious,JPrevious,cfg)

m=numel(z); alpha=cfg.vb.alpha;
jPrior=cfg.vb.rho*(jPrevious-m-1)+m+1;
JPrior=cfg.vb.rho*JPrevious;
x=xp; P=Pp;
for iteration=1:cfg.vb.iterations
    M=(x-xp)*(x-xp)'+P;
    Peff=regularize_cov((alpha*Pp+M)/(alpha+1),cfg.epsilon);
    residual=z-H*x;
    O=residual*residual'+H*P*H';
    j=jPrior+1; J=JPrior+O;
    R=regularize_cov(J/(j-m-1),cfg.epsilon);
    [x,P]=gaussian_update(xp,Peff,z,H,R,cfg.epsilon);
end

end
