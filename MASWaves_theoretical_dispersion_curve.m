
function [c_t,lambda_t] = MASWaves_theoretical_dispersion_curve(c_test,lambda,h,alpha,beta,rho,n)

% Compute wave numbers that correspond to wavelengths lambda
k = (2*pi)./lambda;

D = zeros(length(c_test),length(k));
c_t = zeros(length(k),1);

% For each wave number, recompute the system stiffness matrix using
% different values of c_test until its determinant has a sign change. 
for l = 1:length(k)
    for m = 1:length(c_test)
        D(l,m) = MASWaves_stiffness_matrix(c_test(m),k(l),h,alpha,beta,rho,n);
        if m==1;
            sign_old = sign(D(l,m));
        else
            sign_old = signD;
        end
        signD = sign(D(l,m));
        if sign_old*signD == -1
            c_t(l)=c_test(m);
            lambda_t(l)=2*pi/k(l);
            break
        end
    end
end
end