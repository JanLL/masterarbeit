% rational splines testing

nrb_order = 4; % nrb_order = 4 equates to C^2

cntrl_pts = [0, 30, 60, 90, 120, 125, 130, 132., 135, 150, 160, 180; ...
             1., 1,  1.1, 1.15, 1.2, 5., 10, 1.5, 1.51, 1.52, 1.53, 1.54];
% cntrl_pts = [30, 60, 90, 120, 125, 130, 130., 135, 150, 160; ...
%              0, 0, 0, 0, 0, 1, 0, 0, 0, 0];         
         
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
C = nrbeval(curve, tt);

c_p = interp1(C(1,:), C(2,:), 'linear', 'pp');
dc_p = fnder(c_p, 1);

T = 30:0.01:160;
plot(T, ppval(c_p, T)); hold on
plot(T, ppval(dc_p, T)); 



return;




tt = 0:0.001:1;
[C,dC] = nrbdeval(curve, dcurve, tt);

x = C(1,:);
y = C(2,:);


dy = dC(2,:) ./ dC(1,:); % gradients -> derivative

x_interp = 30:0.005:160;

tic;
y_interp = interp1(x, y, x_interp);
dy_interp = interp1(x, dy, x_interp);
toc

figure();
%cla;

plot(x_interp, y_interp, 'DisplayName', 'Matlab'); hold on

%plot(x_interp, dy_interp)


data_nurbes = importdata('/home/argo/masterarbeit/code/tools/nurbs_c++/curve_nurbes.txt');
plot(data_nurbes(:,1), data_nurbes(:,2), '--', 'DisplayName', 'NURBES')

data_spline = importdata('/home/argo/masterarbeit/code/tools/nurbs_c++/curve_interpol.txt');
plot(data_spline(:,1), data_spline(:,2), '--', 'DisplayName', 'Spline')

legend('show', 'location', 'northwest')

