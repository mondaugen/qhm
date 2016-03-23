function [x] = synth_linear_amod_sinusoid(a0,a1,f0,phi0=0,Fs=8000,T=1.0)
t=(0:(1/Fs):T);
m=(a1-a0)/T;
x=(m*t+a0).*exp(1i*(2*pi*f0*t+phi0));
