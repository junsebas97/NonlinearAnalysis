function [C, S_vec, S_mat] = hyperelast(lambda, mu, F)
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
ln_f = log(det(F));

% the coefficients of the constitutive matrix are calculated and it is
% assembled
C_11 = (-2*ln_f*lambda + lambda + 2*mu)*c(1, 1)^2;
C_12 = (2*mu - 2*ln_f*lambda)*c(1, 2)^2 + lambda*c(1, 1)*c(2, 2);
C_13 = (2*mu - 2*ln_f*lambda)*c(3, 1)^2 + lambda*c(1, 1)*c(3, 3);
C_14 = (-2*ln_f*lambda + lambda + 2*mu)*c(1, 2)*c(1, 1);
C_15 = lambda*c(1, 1)*c(2, 3) + 2*(mu - ln_f*lambda)*c(1, 2)*c(3, 1);
C_16 = (-2*ln_f*lambda + lambda + 2*mu)*c(1, 1)*c(3, 1);
C_22 = (-2*ln_f*lambda + lambda + 2*mu)*c(2, 2)^2;
C_23 = (2*mu - 2*ln_f*lambda)*c(2, 3)^2 + lambda*c(2, 2)*c(3, 3);
C_24 = (-2*ln_f*lambda + lambda + 2*mu)*c(1, 2)*c(2, 2);
C_25 = (-2*ln_f*lambda + lambda + 2*mu)*c(2, 2)*c(2, 3);
C_26 = (2*mu - 2*ln_f*lambda)*c(1, 2)*c(2, 3) + lambda*c(2, 2)*c(3, 1);
C_33 = (-2*ln_f*lambda + lambda + 2*mu)*c(3, 3)^2;
C_34 = (2*mu - 2*ln_f*lambda)*c(2, 3)*c(3, 1) + lambda*c(1, 2)*c(3, 3);
C_35 = (-2*ln_f*lambda + lambda + 2*mu)*c(2, 3)*c(3, 3);
C_36 = (-2*ln_f*lambda + lambda + 2*mu)*c(3, 1)*c(3, 3);
C_44 = (-ln_f*lambda + lambda + mu)*c(1, 2)^2 + (mu - ln_f*lambda)*c(1, 1)*c(2, 2);
C_45 = (-ln_f*lambda + lambda + mu)*c(1, 2)*c(2, 3) + (mu - ln_f*lambda)*c(2, 2)*c(3, 1);
C_46 = (-ln_f*lambda + lambda + mu)*c(1, 2)*c(3, 1) + (mu - ln_f*lambda)*c(1, 1)*c(2, 3);
C_55 = (-ln_f*lambda + lambda + mu)*c(2, 3)^2 + (mu - ln_f*lambda)*c(2, 2)*c(3, 3);
C_56 = (-ln_f*lambda + lambda + mu)*c(2, 3)*c(3, 1) + (mu - ln_f*lambda)*c(1, 2)*c(3, 3);
C_66 = (-ln_f*lambda + lambda + mu)*c(3, 1)^2 + (mu - ln_f*lambda)*c(1, 1)*c(3, 3);

C = [C_11, C_12, C_13, C_14, C_15, C_16;
     C_12, C_22, C_23, C_24, C_25, C_26;
     C_13, C_23, C_33, C_34, C_35, C_36;
     C_14, C_24, C_34, C_44, C_45, C_46;
     C_15, C_25, C_35, C_45, C_55, C_56;
     C_16, C_26, C_36, C_46, C_56, C_66];

% the second PK tensor is computed and assembled as a matrix and a vector
S = lambda*ln_f*c + mu*(eye(3) - c);

S_mat = zeros(9);
S_mat(1:3, 1:3) = S;       S_mat(4:6, 4:6) = S;       S_mat(7:9, 7:9) = S;

S_vec = [S(1, 1); S(2, 2); S(3, 3); S(1, 2); S(2, 3); S(1, 3)];
end
