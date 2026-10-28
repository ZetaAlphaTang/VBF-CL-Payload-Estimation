function [C,lambda_reg]=regularize_cov(C,epsilon)
C=(C+C')/2;
lambda_reg=0;
[~,flag]=chol(C);
if flag~=0
    minimum_eigenvalue=min(real(eig(C)));
    lambda_reg=epsilon-minimum_eigenvalue;
    C=C+lambda_reg*eye(size(C,1));
    C=(C+C')/2;
end
end
