q = .3; 
v = 0;
z = 0;

K = 100;
D = 10;

h= .001;
tf = 25;

N= tf/h;

ts = zeros(1, N+1);
qs = zeros(1, N+1);
zs = zeros(1, N+1);

ts(1) = 0;
qs(1) = q;
zs(1) = z;

for i=1:N
  if (q > z)
    f = 0;
  else
    f = max(0, -K*z - D*v);
  end
  dz = - (K*z + f)/D;
  
  if (ts(i)<tf/2)
    u = 0;
  else
    u = 1.1
  end
  
  v = v + h*(f - 1 + u);
  q = q + h*v;
  z = z + h*dz;
  qs(i+1)= q;
  zs(i+1)= z;
  ts(i+1)= i*h;
end

plot(ts, qs, ts, zs);
legend('q','z')