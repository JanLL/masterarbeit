startup;

i0 = find(t > 69, 1, 'first');
%ie = find(T(:,1) < T(:,2), 1, 'first');
ie = find(t < 74, 1, 'last');

i0 = 10000


dt = t(2) - t(1); % assume: equidistant

dT1 = diff(T(:,2))./dt;
dT2 = diff(T(:,3))./dt;
dT5 = diff(T(:,6))./dt;
dTN = diff(T(:,end))./dt;


ddT1 = diff(dT1)./dt;
ddT2 = diff(dT2)./dt;
ddT5 = diff(dT5)./dt;
ddTN = diff(dTN)./dt;


figure(1);
plot(t(i0:ie), ddT1(i0:ie), 'DisplayName', 'ddT_1'); hold on
plot(t(i0:ie), ddT2(i0:ie), 'DisplayName', 'ddT_2'); hold on
plot(t(i0:ie), ddT5(i0:ie), 'DisplayName', 'ddT_5'); hold on

%plot(t(1:ie), ddTN(1:ie), 'DisplayName', 'ddT_N'); hold on

legend(gca, 'show', 'Location', 'northwest')
xlabel('Time t');

 
figure(2);
plot(t(i0:ie), T(i0:ie,1), 'g--', 'DisplayName', 'T_0'); hold on
plot(t(i0:ie), T(i0:ie,2), 'r--', 'DisplayName', 'T_1'); hold on
plot(t(i0:ie), T(i0:ie,6), 'r--', 'DisplayName', 'T_5'); hold on
%plot(t(1:ie), T(1:ie,end), 'r--', 'DisplayName', 'T_N'); hold on

xlabel('Time t');
ylabel('Temp T');