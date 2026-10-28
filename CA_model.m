function [F,G] = CA_model(dt)
validateattributes(dt,{'double'},{'scalar','positive','finite'});
F=[eye(3),dt*eye(3);zeros(3),eye(3)];
G=[0.5*dt^2*eye(3);dt*eye(3)];
end
