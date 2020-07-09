function J = jacobian(coord, GP)
%{
This function calcules the Jacobien between the passed coordinates and the
interpolation space evaluated in a Gauss point of the element

Args:
coord: nodal coordinates of the FE [x, y]
pg:    Gauss point

Outputs: the Jacobian
%}
s = GP(1);        t = GP(2);        r = GP(3);

n_nod = size(coord, 1);

% derivatives of the interpolation functions:
[~, dN_ds, dN_dt, dN_dr, ~, ~] = parameters(n_nod);

dx0_ds = coord(:, 1)'*dN_ds(s, t, r);
dx0_dt = coord(:, 1)'*dN_dt(s, t, r);
dx0_dr = coord(:, 1)'*dN_dr(s, t, r);

dy0_ds = coord(:, 2)'*dN_ds(s, t, r);
dy0_dt = coord(:, 2)'*dN_dt(s, t, r);
dy0_dr = coord(:, 2)'*dN_dr(s, t, r);

dz0_ds = coord(:, 3)'*dN_ds(s, t, r);
dz0_dt = coord(:, 3)'*dN_dt(s, t, r);
dz0_dr = coord(:, 3)'*dN_dr(s, t, r);

% the Jacobian matrix:
J = [dx0_ds, dx0_dt, dx0_dr;
     dy0_ds, dy0_dt, dy0_dr;
     dz0_ds, dz0_dt, dz0_dr];
end
