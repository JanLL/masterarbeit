% rational splines testing

nrb_order = 4; % nrb_order = 4 equates to C^2

cntrl_pts = [30, 60, 90, 120, 125, 130, 131, 135, 150; ...
             1,  1.1, 1.15, 1.2, 5., 10, 1.5, 1.51, 1.52];
num_cntrl_pts = size(cntrl_pts,2);

% equidistant knots
knots = [zeros(1,nrb_order), ...
         (1:num_cntrl_pts-nrb_order)/(num_cntrl_pts-nrb_order+1), ...
         ones(1,nrb_order)];

% non-equidistant knots
% knots = [0.0, 0.0, 0.0, 0.0, 0.2, 0.4, 0.7, 0.8 1.0, 1.0, 1.0, 1.0];


length(knots)
size(cntrl_pts, 2)

%assert(length(knots) - size(cntrl_pts, 2) == 3);

curve = nrbmak(cntrl_pts, knots);
dcurve = nrbderiv(curve);

           
tt = 0:0.01:1;

tic;
[C,dC] = nrbdeval(curve, dcurve, tt);
toc


dx = dC(1,:);
dy = dC(2,:) ./ dC(1,:);


figure(1);
%cla;
plot(x,y); hold on
plot(x, dy)