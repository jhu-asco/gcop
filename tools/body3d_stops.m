function f = body3d_stops

vs = 4:4:20;
ds = [0.931897 2.22661 4.23726 7.41835 12.7978]
fs = [21.5183 28.5696 32.1243 32.1491 30.3235]
ts = [0.4 0.5 0.65 0.85 1.15]

figure
plot(vs, ds, '-.', 'LineWidth',3)
ax = gca;
%ax.XTick = 1:M;
%ax.XTickLabel = names;
title('Stopping Distance vs. Starting Velocity')
xlabel('m/s')
ylabel('m')

figure
plot(vs, ts, '-.', 'LineWidth',3)
ax = gca;
%ax.XTick = 1:M;
%ax.XTickLabel = names;
title('Time to Stop vs Starting Velocity')
xlabel('m/s')
ylabel('s')

figure
plot(vs, fs, '-.', 'LineWidth',3)
ax = gca;
%ax.XTick = 1:M;
%ax.XTickLabel = names;
title('Max Thrust to Stop vs. Starting Velocity')
xlabel('m/s')
ylabel('N')


