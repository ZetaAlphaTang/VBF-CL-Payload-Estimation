function [z,H]=linearized_measurement(xp,obs,k,model)
if strcmp(model,'CL-pseudolinear')
    z=obs.pUav(:,k)+obs.length(k)*obs.bearing(:,k);
    H=[eye(3),zeros(3)]; return;
end
rvec=xp(1:3)-obs.pUav(:,k); r=norm(rvec);
if r<1e-9, error('Undefined bearing: predicted relative range is zero.'); end
b=rvec/r; E=obs.tangent(:,:,k);
H=[E'*(eye(3)-b*b')/r,zeros(2,3)];
z=E'*obs.bearing(:,k)-E'*b+H*xp;

end
