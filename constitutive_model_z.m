function [sigma_extra] = constitutive_model_z(F)

%=========================================================================
% This function calculates the Cauchy stress extra (sigma_extra) for a  
% given gradient tensor F, Using Neo-Hookean Constitutive model (see Lecture 8)
%=========================================================================
global inputDataExp

% Extract material constants from input data and store in vector
c = inputDataExp.estimated_parameters; % material constant, kPa

% Calculate the Right Cauchy-Green Deformation tensor
C = F'*F;

% Define identity matrixd
I = eye(3);

% Calculate the Green-Lagrange Strain tensor
E = 0.5*(C-I);

% Calculate the derivative of W wrt Ezz
Q = (c(1)*E(1,1)^2)+(c(2)*E(2,2)^2)+(c(3)*E(3,3)^2)+(2*c(4)*E(1,1)*E(2,2))...
    +(2*c(5)*E(2,2)*E(3,3))+(2*c(6)*E(1,1)*E(3,3));
dW = (2*c(3)*E(3,3)+2*c(5)*E(2,2)+2*c(6)*E(1,1)).*0.5*c(7)*exp(Q);
sigma_extra = C.*dW; % Cauchy extra stress

end