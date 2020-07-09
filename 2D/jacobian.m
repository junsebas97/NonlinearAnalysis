function J = jacobian(coord, GP)
%{
This function calcules the Jacobien between the passed coordinates and the
interpolation space evaluated in a Gauss point of the element

Args:
coord: nodal coordinates of the FE [x, y]
pg:    Gauss point

Outputs: the Jacobian
%}
s = GP(1);        t = GP(2);

n_nod = size(coord, 1);

% derivatives of the interpolation functions:
if n_nod == 4
    dN_ds = [(t - 1)/4;  (1 - t)/4; (t + 1)/4; -(t + 1)/4];
    dN_dt = [(s - 1)/4; -(s + 1)/4; (s + 1)/4;  (1 - s)/4];
elseif n_nod == 8
    dN_ds = [-((2*s + t)*(t - 1))/4;           s*(t - 1);
             -((2*s - t)*(t - 1))/4;         1/2 - t^2/2;
              ((2*s + t)*(t + 1))/4;          -s*(t + 1);
              ((2*s - t)*(t + 1))/4;         t^2/2 - 1/2];

    dN_dt = [-((s + 2*t)*(s - 1))/4;         s^2/2 - 1/2;
             -((s - 2*t)*(s + 1))/4;          -t*(s + 1);
              ((s + 2*t)*(s + 1))/4;         1/2 - s^2/2;
              ((s - 2*t)*(s - 1))/4;           t*(s - 1)];
end

dx0_ds = coord(:, 1)'*dN_ds;        dx0_dt = coord(:, 1)'*dN_dt;
dy0_ds = coord(:, 2)'*dN_ds;        dy0_dt = coord(:, 2)'*dN_dt;

% the Jacobian matrix:
J = [dx0_ds, dx0_dt; dy0_ds, dy0_dt];
end
