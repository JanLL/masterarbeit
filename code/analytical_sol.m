function [T] = analytical_sol(x, t, n, T0, heat_rate, a)
% Analytical solution T(x,t) of heat equation with boundary conditions: 
% T(x(1)  ,t)            = T0 + heat_rate*t
% \partial_x T(x(end),t) = 0
% Initial condition: T(x,0) = T0
%
% INPUT:
%         x --> spatial domain of solution.
%         t --> time domain of solution.
%         n --> number of (actually infinite) sum elements.
%        T0 --> Initial temperature: T(x,0) = T0
% heat_rate --> heat rate [K/s] of Dirichlet boundary condition 
%               T(x(1)  ,t = T0 + heat_rate*t
%         a --> temperature conductivity
%
% OUTPUT:
%       T --> 2D array with dimensions [length(t), length(x)] with
%             temperature at point (x,t)

N = 1:n;
N = reshape(N,[1,1,length(N)]);
N = repmat(N,[length(t), length(x), 1]);

[x_grid, t_grid] = meshgrid(x,t);
x_grid = repmat(x_grid, [1,1,size(N,3)]);
t_grid = repmat(t_grid, [1,1,size(N,3)]);

L = x(end);
lambda_sqrt = @(n) ((2*n-1)*pi/(2*L)); % optimierungspotential nur einmal lambda auszurechnen...



T = T0 + heat_rate*t_grid(:,:,1) - sum(4.*heat_rate./(pi.*(2.*N-1)) .* 1./(a.*lambda_sqrt(N).^2) .* ...
          (1 - exp(-a.*lambda_sqrt(N).^2.*t_grid)) .* sin(lambda_sqrt(N).*x_grid),3);

      
      
      
end

