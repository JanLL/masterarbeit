function [u] = analytical_sol(x, t, n, T0, heat_rate, a)
%

N = 1:n;
N = reshape(N,[1,1,length(N)]);
N = repmat(N,[length(t), length(x), 1]);

[X, T] = meshgrid(x,t);
X = repmat(X, [1,1,size(N,3)]);
T = repmat(T, [1,1,size(N,3)]);

L = x(end);
lambda_sqrt = @(n) ((2*n-1)*pi/(2*L)); % optimierungspotential nur einmal lambda auszurechnen...



u = T0 + heat_rate*T(:,:,1) - sum(4.*heat_rate./(pi.*(2.*N-1)) .* 1./(a.*lambda_sqrt(N).^2) .* ...
          (1 - exp(-a.*lambda_sqrt(N).^2.*T)) .* sin(lambda_sqrt(N).*X),3);

      
      
      
end

