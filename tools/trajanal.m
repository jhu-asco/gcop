function f = trajanal

%type 1:RJET, 2:MPPT, 3:VACCO

mass = 1

if mass

  names = {'1.25 kg', '1.5 kg', '1.75 kg', '2 kg'};
  fnames = {'../logs/body3d/runs/mass_env/traj4_1.txt',...
            '../logs/body3d/runs/mass_env/traj4_2.txt',...
            '../logs/body3d/runs/mass_env/traj4_3.txt',...
            '../logs/body3d/runs/mass_env/traj4_4.txt'};

else

names = {'Empty', 'Slalom1 (easy)', 'Slalom2 (mod.)', 'Slalom3 (hard)'};
fnames = {'../logs/body3d/runs/empty/traj.txt', ...
          '../logs/body3d/runs/slalom1/traj.txt',...
          '../logs/body3d/runs/slalom2/traj.txt',...
          '../logs/body3d/runs/slalom3/traj.txt'};

end

M = length(names)

vmax = zeros(M,1);
vave = zeros(M,1);
fmax = zeros(M,1);
fave = zeros(M,1);

ttot = zeros(M,1);

W = zeros(M,1);

res = 5;

for i=1:M
D = load(fnames{i}, 'ascii')';

ts = D(1,1:res:end-1);
vs = D(2:7,1:res:end-1);
qs = D(8:13,1:res:end-1);
us = D(14:end,1:res:end-1);

W(i) = 0;
for j=1:length(ts)-1;
  dt = ts(j+1)-ts(j);
  W(i) = W(i) + dt*(norm(vs(4:6,j))*us(4,j) + vs(1:3,i)'*us(1:3,i));
end
Pave = W(i)/(ts(end)-ts(1));

dt = ts(2)-ts(1);

vns = sqrt(sum(vs(4:6,:).*vs(4:6,:), 1))';  % transl vel norms

ttot(i) = ts(end) - ts(1);
vmax(i) = max(vns);
vave(i) = mean(vns);

fmax(i) = max(us(4,:));
fave(i) = mean(us(4,:));

fws = sqrt(sum(us(1:3,:).*us(1:3,:), 1))';  % transl vel norms
fwmax(i) = max(fws);
fwave(i) = mean(fws);


end

figure
plot(1:M, vave, 1:M, vmax, '-.', 'LineWidth',3)
ax = gca;
ax.XTick = 1:M;
ax.XTickLabel = names;
legend('V_{ave}', 'V_{max}')
ylabel('m/s')

figure
plot(1:M, fave, 1:M, fmax, '-.', 'LineWidth',3)
ax = gca;
ax.XTick = 1:M;
ax.XTickLabel = names;
legend('F_{ave}', 'F_{max}')
ylabel('N')

figure
plot(1:M, fwave, 1:M, fwmax, '-.', 'LineWidth',3)
ax = gca;
ax.XTick = 1:M;
ax.XTickLabel = names;
legend('|\tau|_{ave}', '|\tau|_{max}')
ylabel('N/m')


figure
plot(1:M, ttot, 'LineWidth',3)
ax = gca;
ax.XTick = 1:M;
ax.XTickLabel = names;
legend('T_{total}')
ylabel('sec')


figure
plot(1:M, W, 'LineWidth',3)
ax = gca;
ax.XTick = 1:M;
ax.XTickLabel = names;
legend('W_{total}')
ylabel('sec')
