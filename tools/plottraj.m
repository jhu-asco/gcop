function f = plottraj(fname, type)

%type 1:RJET, 2:MPPT, 3:VACCO

res = 5;

D = load(fname, 'ascii')';

ts = D(1,1:res:end-1)
vs = D(2:7,1:res:end-1);
qs = D(8:13,1:res:end-1);
us = D(14:end,1:res:end-1);

W = 0;
for i=1:length(ts)-1;
  dt = ts(i+1)-ts(i);
  W = W + dt*norm(vs(4:6,i))*us(4,i);
%  W = W + dt*vs(1:3,i)'*us(1:3,i);
end
Pave = W/(ts(end)-ts(1))

display(['Motor Work= ' num2str(W) ' Watts-seconds']);
display(['Motor Power (ave)= ' num2str(Pave) ' Watts']);


%return

figure
plot(ts, vs(1,:), '-r', ts, vs(2,:), '-og', ts, vs(3,:), '-db', 'LineWidth',3);

set(gca, 'FontSize',25)
xlabel('sec.')
ylabel('rad/s')
title('Angular velocity')
legend( '\omega_x','\omega_y','\omega_z')

%saveas(gca,['w.eps'],'psc2');

figure
plot(ts, vs(4,:), '-r', ts, vs(5,:), '-og', ts,  vs(6,:), '-db', 'LineWidth',3)

legend( 'v_x','v_y','v_z')
set(gca, 'FontSize',25)
xlabel('sec.')
ylabel('m/s')
title('Linear velocity')
%saveas(gca,['v.eps'],'psc2');

dt = ts(2)-ts(1);

vns = sqrt(sum(vs(4:6,:).*vs(4:6,:), 1));  % transl vel norms
vmax = max(vns);
vmean = mean(vns);
hold on

plot(ts, vmax*ones(size(ts)), ts, vmean*ones(size(ts)), 'LineWidth',4)
legend( 'v_x','v_y','v_z','v_{max}','v_{mean}')

figure

h = plot(ts, us','LineWidth',3)

if(size(us,1)==6)
  set(h(1),'LineStyle','--');
  set(h(2),'LineStyle','-.');
  set(h(3),'Marker','+');
  set(h(4),'Marker','o');
  set(h(5),'LineStyle','-');
  set(h(5),'LineStyle','-o');
end

if(size(us,1)==4)
  set(h(1),'LineStyle','--');
  set(h(2),'LineStyle','-.');
  set(h(3),'Marker','+');
  set(h(4),'LineStyle','-');
end


disp('total delta V (m/s)')  
sum(sum(us*dt))/2

c = size(us, 1);
switch (c)
 case 4
  legend('u_1','u_2','u_3','u_4')
 case 6
  legend('u_1','u_2','u_3','u_4','u_5','u_6')
  
 otherwise
end

set(gca, 'FontSize',25)
title('Control inputs');

%legend('u_\omega','u_r')
xlabel('sec.')
ylabel('N or N/m')
%saveas(gca,['u.eps'],'psc2');
