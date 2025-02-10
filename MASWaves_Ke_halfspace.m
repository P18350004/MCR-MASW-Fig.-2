
function Ke_halfspace = MASWaves_Ke_halfspace(alpha,beta,rho,c_test,k)

r = sqrt(1-c_test^2/alpha^2);
s = sqrt(1-c_test^2/beta^2);

k_11 = k*rho*beta^2*(r*(1-s^2))/(1-r*s);
k_12 = k*rho*beta^2*(1-s^2)/(1-r*s) - 2*k*rho*beta^2;
k_21 = k_12;
k_22 = k*rho*beta^2*(s*(1-s^2))/(1-r*s);

Ke_halfspace = [k_11 k_12 ; k_21 k_22];

end