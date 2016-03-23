T=1.0;
Fs=8000;
f0=350;
phi0=0.5;
t=(0:1/Fs:T);
y=exp(1j*(2*pi*f0*t + phi0));%+0.001*rand(1,length(t));
y=y(:);
b=[0.5; 350; 1;];
U=1;
for u=(1:U);
    g=exp(1j*(2*pi*b(2)*t'+b(1)));
    J=[-1j*b(3)*g -1j*2*pi*t'*b(3).*g -g];
    cond(J'*J)
    b=b+(J'*J)^(-1)*(J'*(y-b(3)*g));
end
plot(t,y,'.',t,b(3)*exp(1j*(2*pi*b(2)*t'+b(1))),'-');
abs(b)
