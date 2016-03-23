T=1.0;
Fs=8000;
f0=350;
phi0=0.5;
t=(0:1/Fs:T);
y=exp(1j*(2*pi*f0*t + phi0));%+0.001*rand(1,length(t));
y=y(:);
b=[0.5; 350; 1;];
U=100;
gamma=0.1;
for u=(1:U);
    g=1j*(2*pi*b(2)*t'+b(1));
    delta=[sum(-2*b(3)*exp(g).*(y-b(3)*exp(g))); ...
           sum(-2*b(3)*1j*2*pi*t'.*exp(g).*(y-b(3)*exp(g))); ...
           sum(-2*exp(g).*(y-b(3)*exp(g)))];
    b=b-gamma*delta;
end
plot(t,y,'.',t,b(3)*exp(1j*(2*pi*b(2)*t'+b(1))),'-');
b
