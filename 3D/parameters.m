function [N, dN_ds, dN_dt, dN_dr, pg, w] = parameters(number_nodes)
%{
This function returns the interpolation functions and the Gauss quadrature
of the FE

Args: Number of nodes of the FE

Returns:
N:      interpolation functions
dN_ds:  derivative of the interpolation function w.r.t to S-axis
dN_dt:  derivative of the interpolation function w.r.t to T-axis
dN_dr:  derivative of the interpolation function w.r.t to R-axis
pg:     coordinates of the Gauss points 
weight: weights of the Gauss points
%}
if number_nodes == 20
    % interpolation functions and their derivatives:
    N = @(s, t, r)[((t - 1)*(s - 1)*(r - 1)*(t + s + r + 2))/8;     % N1
                                -((s^2 - 1)*(t - 1)*(r - 1))/4;     % N2
                  -((t - 1)*(s + 1)*(r - 1)*(t - s + r + 2))/8;     % N3
                                 ((t^2 - 1)*(s + 1)*(r - 1))/4;     % N4
                  -((t + 1)*(s + 1)*(r - 1)*(t + s - r - 2))/8;     % N5
                                 ((s^2 - 1)*(t + 1)*(r - 1))/4;     % N6
                  -((t + 1)*(s - 1)*(r - 1)*(s - t + r + 2))/8;     % N7
                                -((t^2 - 1)*(s - 1)*(r - 1))/4;     % N8
                                -((r^2 - 1)*(t - 1)*(s - 1))/4;     % N9
                                 ((r^2 - 1)*(t - 1)*(s + 1))/4;     % N10
                                -((r^2 - 1)*(t + 1)*(s + 1))/4;     % N11
                                 ((r^2 - 1)*(t + 1)*(s - 1))/4;     % N12
                  -((t - 1)*(s - 1)*(r + 1)*(t + s - r + 2))/8;     % N13
                                 ((s^2 - 1)*(t - 1)*(r + 1))/4;     % N14
                   ((t - 1)*(s + 1)*(r + 1)*(t - s - r + 2))/8;     % N15
                                -((t^2 - 1)*(s + 1)*(r + 1))/4;     % N16
                   ((t + 1)*(s + 1)*(r + 1)*(t + s + r - 2))/8;     % N17
                                -((s^2 - 1)*(t + 1)*(r + 1))/4;     % N18
                  -((t + 1)*(s - 1)*(r + 1)*(t - s + r - 2))/8;     % N19
                                 ((t^2 - 1)*(s - 1)*(r + 1))/4];    % N20

    dN_ds = @(s, t, r)[((r - 1)*(t - 1)*(r + 2*s + t + 1))/8;     % dN1/ds
                                      -(s*(r - 1)*(t - 1))/2;     % dN2/ds
                      -((r - 1)*(t - 1)*(r - 2*s + t + 1))/8;     % dN3/ds
                                       ((t^2 - 1)*(r - 1))/4;     % dN4/ds
                       ((r - 1)*(t + 1)*(r - 2*s - t + 1))/8;     % dN5/ds
                                       (s*(r - 1)*(t + 1))/2;     % dN6/ds
                      -((r - 1)*(t + 1)*(r + 2*s - t + 1))/8;     % dN7/ds
                                      -((t^2 - 1)*(r - 1))/4;     % dN8/ds
                                      -((r^2 - 1)*(t - 1))/4;     % dN9/ds
                                       ((r^2 - 1)*(t - 1))/4;     % dN10/ds
                                      -((r^2 - 1)*(t + 1))/4;     % dN11/ds
                                       ((r^2 - 1)*(t + 1))/4;     % dN12/ds
                      -((r + 1)*(t - 1)*(2*s - r + t + 1))/8;     % dN13/ds
                                       (s*(r + 1)*(t - 1))/2;     % dN14/ds
                      -((r + 1)*(t - 1)*(r + 2*s - t - 1))/8;     % dN15/ds
                                      -((t^2 - 1)*(r + 1))/4;     % dN16/ds
                       ((r + 1)*(t + 1)*(r + 2*s + t - 1))/8;     % dN17/ds
                                      -(s*(r + 1)*(t + 1))/2;     % dN18/ds  
                      -((r + 1)*(t + 1)*(r - 2*s + t - 1))/8;     % dN19/ds
                                       ((t^2 - 1)*(r + 1))/4];    % dN20/ds

    dN_dt = @(s, t, r)[ ((r - 1)*(s - 1)*(r + s + 2*t + 1))/8;    % dN1/dt
                                       -((s^2 - 1)*(r - 1))/4;    % dN2/dt
                       -((r - 1)*(s + 1)*(r - s + 2*t + 1))/8;    % dN3/dt
                                        (t*(r - 1)*(s + 1))/2;    % dN4/dt
                        ((r - 1)*(s + 1)*(r - s - 2*t + 1))/8;    % dN5/dt
                                        ((s^2 - 1)*(r - 1))/4;    % dN6/dt
                       -((r - 1)*(s - 1)*(r + s - 2*t + 1))/8;    % dN7/dt
                                       -(t*(r - 1)*(s - 1))/2;    % dN8/dt
                                       -((r^2 - 1)*(s - 1))/4;    % dN9/dt
                                        ((r^2 - 1)*(s + 1))/4;    % dN10/dt
                                       -((r^2 - 1)*(s + 1))/4;    % dN11/dt
                                        ((r^2 - 1)*(s - 1))/4;    % dN12/dt
                       -((r + 1)*(s - 1)*(s - r + 2*t + 1))/8;    % dN13/dt
                                        ((s^2 - 1)*(r + 1))/4;    % dN14/dt
                       -((r + 1)*(s + 1)*(r + s - 2*t - 1))/8;    % dN15/dt
                                       -(t*(r + 1)*(s + 1))/2;    % dN16/dt
                        ((r + 1)*(s + 1)*(r + s + 2*t - 1))/8;    % dN17/dt
                                       -((s^2 - 1)*(r + 1))/4;    % dN18/dt
                       -((r + 1)*(s - 1)*(r - s + 2*t - 1))/8;    % dN19/dt
                                        (t*(r + 1)*(s - 1))/2];   % dN20/dt
                                    
	dN_dr = @(s, t, r)[ ((s - 1)*(t - 1)*(2*r + s + t + 1))/8;     %dN1/dr
                                       -((s^2 - 1)*(t - 1))/4;     %dN2/dr
                       -((s + 1)*(t - 1)*(2*r - s + t + 1))/8;     %dN3/dr
                                        ((t^2 - 1)*(s + 1))/4;     %dN4/dr
                        ((s + 1)*(t + 1)*(2*r - s - t + 1))/8;     %dN5/dr
                                        ((s^2 - 1)*(t + 1))/4;     %dN6/dr
                       -((s - 1)*(t + 1)*(2*r + s - t + 1))/8;     %dN7/dr
                                       -((t^2 - 1)*(s - 1))/4;     %dN8/dr
                                       -(r*(s - 1)*(t - 1))/2;     %dN9/dr
                                        (r*(s + 1)*(t - 1))/2;     %dN10/dr
                                       -(r*(s + 1)*(t + 1))/2;     %dN11/dr
                                        (r*(s - 1)*(t + 1))/2;     %dN12/dr
                       -((s - 1)*(t - 1)*(s - 2*r + t + 1))/8;     %dN13/dr
                                        ((s^2 - 1)*(t - 1))/4;     %dN14/dr
                       -((s + 1)*(t - 1)*(2*r + s - t - 1))/8;     %dN15/dr
                                       -((t^2 - 1)*(s + 1))/4;     %dN16/dr
                        ((s + 1)*(t + 1)*(2*r + s + t - 1))/8;     %dN17/dr
                                       -((s^2 - 1)*(t + 1))/4;     %dN18/dr
                       -((s - 1)*(t + 1)*(2*r - s + t - 1))/8;     %dN19/dr
                                        ((t^2 - 1)*(s - 1))/4];    %dN20/dr

    % Gauss-Legendre quadrature:
    pg = sqrt(1/3)*[-1, -1, -1; 1, -1, -1; 1, 1, -1; -1, 1, -1;...
                    -1, -1,  1; 1, -1,  1; 1, 1,  1; -1, 1,  1];
    w  = [1, 1, 1, 1, 1, 1, 1, 1];
end
end