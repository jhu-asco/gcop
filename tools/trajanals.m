function f = trajanals

%type 1:RJET, 2:MPPT, 3:VACCO

mass = 1

  names = {'1.2 kg', '1.5 kg', '1.8 kg', '2.1 kg'};
  fnames = {{'../logs/body3d/runs/mass_env/traj1_1.txt',...
            '../logs/body3d/runs/mass_env/traj1_2.txt',...
            '../logs/body3d/runs/mass_env/traj1_3.txt',...
            '../logs/body3d/runs/mass_env/traj1_4.txt'},
            {'../logs/body3d/runs/mass_env/traj2_1.txt',...
            '../logs/body3d/runs/mass_env/traj2_2.txt',...
            '../logs/body3d/runs/mass_env/traj2_3.txt',...
            '../logs/body3d/runs/mass_env/traj2_4.txt'},
            {'../logs/body3d/runs/mass_env/traj3_1.txt',...
            '../logs/body3d/runs/mass_env/traj3_2.txt',...
            '../logs/body3d/runs/mass_env/traj3_3.txt',...
            '../logs/body3d/runs/mass_env/traj3_4.txt'},
            {'../logs/body3d/runs/mass_env/traj4_1.txt',...
            '../logs/body3d/runs/mass_env/traj4_2.txt',...
            '../logs/body3d/runs/mass_env/traj4_3.txt',...
            '../logs/body3d/runs/mass_env/traj4_4.txt'}};


f1 = figure
plot(0,0)
drawnow
f2 = figure
drawnow
f3 = figure
drawnow
f4 = figure
drawnow
f5 = figure
drawnow

for k=1:4

  M = length(names)
  
  vmax = zeros(M,1);
  vave = zeros(M,1);
  fmax = zeros(M,1);
  fave = zeros(M,1);
  
  ttot = zeros(M,1);
  
  W = zeros(M,1);
  
  res = 5;
  
  for i=1:M
    D = load(fnames{k}{i}, 'ascii')';
    
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
  
figure(f1)
plot(1:M, vave, '-.', 'LineWidth',3)
ax = gca;
ax.XTick = 1:M;
ax.XTickLabel = names;
%legend('V_{ave}')
ylabel('m/s')
hold on

figure(f2)
plot(1:M, fave, '-.', 'LineWidth',3)
ax = gca;
ax.XTick = 1:M;
ax.XTickLabel = names;
% legend('F_{ave}', 'F_{max}')
ylabel('N')
hold on

figure(f3)
plot(1:M, fwave, '-.', 'LineWidth',3)
ax = gca;
ax.XTick = 1:M;
ax.XTickLabel = names;
%legend('|\tau|_{ave}', '|\tau|_{max}')
ylabel('N/m')
hold on


figure(f4)
plot(1:M, ttot, 'LineWidth',3)
ax = gca;
ax.XTick = 1:M;
ax.XTickLabel = names;
legend('T_{total}')
ylabel('sec')
hold on


figure(f5)
plot(1:M, W, 'LineWidth',3)
ax = gca;
ax.XTick = 1:M;
ax.XTickLabel = names;
legend('W_{total}')
ylabel('sec')
hold on

end