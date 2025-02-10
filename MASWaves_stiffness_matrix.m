
function D = MASWaves_stiffness_matrix(c_test,k,h,alpha,beta,rho,n)

% System stiffness matrix
K = zeros(2*(n+1),2*(n+1));

% Check to see if the trial phase velocity is equal to the shear wave velocity
% or compression wave velocity of one of the layers
epsilon = 0.0001;
while any(abs(c_test-beta)<epsilon) || any(abs(c_test-alpha)<epsilon)
    c_test = c_test*(1-epsilon);
end

% Finite thickness layers j = 1,...,n
for j = 1:n
    % Compute element stiffness matrix for layer j
    Ke = MASWaves_Ke_layer(h(j),alpha(j),beta(j),rho(j),c_test,k);
    % Add to the system stiffness matrix
    DOFS = [(2*j-1):(2*j+2)];
    K(DOFS,DOFS) = K(DOFS,DOFS)+Ke;
end

% Half space
% Compute element stiffness matrix for the half space
Ke_halfspace = MASWaves_Ke_halfspace(alpha(end),beta(end),rho(end),c_test,k);
% Add to the system stiffness matrix
DOFS = [2*(n+1)-1:2*(n+1)];
K(DOFS,DOFS) = K(DOFS,DOFS)+Ke_halfspace;

% Evaluate determinant of system stiffness matrix
D = real(det(K));
end