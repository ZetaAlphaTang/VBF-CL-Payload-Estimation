function [x,P]=gaussian_update(xp,Pp,z,H,R,epsilon)
S=regularize_cov(H*Pp*H'+R,epsilon);
K=(Pp*H')/S; x=xp+K*(z-H*xp);
A=eye(6)-K*H; P=regularize_cov(A*Pp*A'+K*R*K',epsilon);
end
