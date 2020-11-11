clc;
clear;
nt = 10;
phi0 = pi/2;
phi = 3*pi/4;
Omega = cos(phi0);
Dt = 1/2;

theta = 0:0.01:2*pi;
rho = sin(theta) - cos(theta);
figure;
polar(theta,rho)
%plot(theta,rho)

e1=[];
e2=[];
for a=1:nt
   tmp1 = exp(-1i*2*pi*(a-1)*Dt*Omega);
   tmp2 = exp(-1i*2*pi*(a-1)*Dt*cos(phi));
   e1 = [e1;    tmp1];
   e2 = [e2;    tmp2];
end
er0 = 1/sqrt(nt) .* e1;
er = 1/sqrt(nt) .* e2;

rho = er0' * er;

figure
polar(abs(rho),'r-o')







