function [N, dN_ds, dN_dt, pg, w] = parameters(number_nodes)
%{
This function returns the interpolation functions and the Gauss quadrature
for the FE with a given number of nodes

Args: Number of nodes of the FE

Returns:
N:      interpolation functions
dN_ds:  derivative of the interpolation function w.r.t to S-axis
dN_ds:  derivative of the interpolation function w.r.t to t-axis
pg:     coordinates of the Gauss points 
weight: weights of the Gauss points
%}
if number_nodes == 4
    % interpolation functions and their derivatives:
    N     = @(s, t)[ ((s - 1)*(t - 1))/4;     % N1
                    -((s + 1)*(t - 1))/4;     % N2
                     ((s + 1)*(t + 1))/4;     % N3
                    -((s - 1)*(t + 1))/4];    % N4

    dN_ds = @(s, t)[ (t - 1)/4;               % dN1/ds
                     (1 - t)/4;               % dN2/ds
                     (t + 1)/4;               % dN3/ds
                    -(t + 1)/4];              % dN4/ds

    dN_dt = @(s, t)[ (s - 1)/4;               % dN1/dt
                    -(s + 1)/4;               % dN2/dt
                     (s + 1)/4;               % dN3/dt
                     (1 - s)/4];              % dN4/dt

    % Gauss-Legendre quadrature:
    pg = sqrt(1/3)*[-1, -1; -1, 1; 1, -1; 1, 1];    % Gauss points
    w  = [1, 1, 1, 1];                              % weights

elseif number_nodes == 8
    N = @(s, t)[-((t - 1)*(s - 1)*(t + s + 1))/4;     % N1
                           ((s^2 - 1)*(t - 1))/2;     % N2
                 ((t - 1)*(s + 1)*(t - s + 1))/4;     % N3
                          -((t^2 - 1)*(s + 1))/2;     % N4
                 ((t + 1)*(s + 1)*(t + s - 1))/4;     % N5
                          -((s^2 - 1)*(t + 1))/2;     % N6
                 ((t + 1)*(s - 1)*(s - t + 1))/4;     % N7
                           ((t^2 - 1)*(s - 1))/2];    % N8

    dN_ds = @(s, t)[-((2*s + t)*(t - 1))/4;           %dN1/ds
                                 s*(t - 1);           %dN2/ds
                    -((2*s - t)*(t - 1))/4;           %dN3/ds
                               1/2 - t^2/2;           %dN4/ds
                     ((2*s + t)*(t + 1))/4;           %dN5/ds
                                -s*(t + 1);           %dN6/ds
                     ((2*s - t)*(t + 1))/4;           %dN7/ds
                               t^2/2 - 1/2];          %dN8/ds

    dN_dt = @(s, t)[-((s + 2*t)*(s - 1))/4;           %dN1/dt
                               s^2/2 - 1/2;           %dN2/dt
                    -((s - 2*t)*(s + 1))/4;           %dN3/dt
                                -t*(s + 1);           %dN4/dt
                     ((s + 2*t)*(s + 1))/4;           %dN5/dt
                               1/2 - s^2/2;           %dN6/dt
                     ((s - 2*t)*(s - 1))/4;           %dN7/dt
                                 t*(s - 1)];          %dN8/dt

    % Gauss-Legendre quadrature
    pg = sqrt(1/3)*[-1, -1; 1, -1; 1, 1; -1, 1];    % Gauss points
    w  = [1, 1, 1, 1];                              % weights
end
end

