% rational splines testing

nrb_order = 4; % nrb_order = 4 equates to C^2

cntrl_pts = [30, 60, 90, 120, 125, 130, 130., 135, 150, 160; ...
             1,  1.1, 1.15, 1.2, 5., 10, 1.5, 1.51, 1.52, 1.53];
num_cntrl_pts = size(cntrl_pts,2);

% equidistant knots
knots = [zeros(1,nrb_order), ...
         (1:num_cntrl_pts-nrb_order)/(num_cntrl_pts-nrb_order+1), ...
         ones(1,nrb_order)];

% non-equidistant knots
% knots = [0.0, 0.0, 0.0, 0.0, 0.2, 0.4, 0.7, 0.8 1.0, 1.0, 1.0, 1.0];


assert(length(knots) - size(cntrl_pts, 2) == 4); % check for C^2 continuity

curve = nrbmak(cntrl_pts, knots);
dcurve = nrbderiv(curve);
           
tt = 0:0.001:1;
tic;
[C,dC] = nrbdeval(curve, dcurve, tt);
toc

x = C(1,:);
y = C(2,:);

dy = dC(2,:) ./ dC(1,:); % gradients -> derivative

x_interp = 30:0.005:160;

tic;
y_interp = interp1(x, y, x_interp);
dy_interp = interp1(x, dy, x_interp);
toc

figure(1);
%cla;

plot(x_interp, y_interp); hold on
plot(x_interp, dy_interp)
