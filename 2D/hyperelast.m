function [C, S_vec, S_mat] = hyperelast(lambda, mu, gamma, F)
%{
Function for Hyperelasticity constitutive matrix and second PK tensor

Args:
lambda, mu: Lamme's constant
gamma:      hyperelasticity constant
F:          deformation gradient

Returns:
C: constitutive matrix for hyperelastic materials
S_vec: equivalent vector of second PK tensor
S_mat: equivalent matrix of second PK tensor
%}

c    = inv(F'*F);            % inverse of the green-cauchy tensor
ln_f = log(det(F)^gamma);

% the coefficients of the constitutive matrix are calculated and it is
% assembled
C_11 = (-2*ln_f*lambda + lambda + 2*mu)*c(1, 1)^2;
C_12 = (2*mu - 2*ln_f*lambda)*c(1, 2)^2 + lambda*c(1, 1)*c(2, 2);
C_13 = (-2*ln_f*lambda + lambda + 2*mu)*c(1, 1)*c(1, 2);
C_22 = (-2*ln_f*lambda + lambda + 2*mu)*c(2, 2)^2;
C_23 = (-2*ln_f*lambda + lambda + 2*mu)*c(1, 2)*c(2, 2);
C_33 = (-ln_f*lambda   + lambda + mu)*c(1, 2)^2 + (mu - ln_f*lambda)*c(1, 1)*c(2, 2);


C = [C_11, C_12, C_13;
     C_12, C_22, C_23;
     C_13, C_23, C_33];

% the second PK tensor is computed and assembled as a matrix and a vector
S = lambda*ln_f*c + mu*(eye(2) - c);
S_mat = zeros(4);         S_mat(1:2, 1:2) = S;         S_mat(3:4, 3:4) = S;
S_vec = [S(1, 1); S(2, 2); S(1, 2)];
end
