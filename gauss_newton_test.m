K=100;
x=randn(K,1);
y=exp(1.8*x)+randn(K,1);
%plot(x,y,'.');
b=1;
U=100;
for u=(1:U),
    J=-b*exp(b*x);
    b=b-(J'*J)^(-1)*J'*(y-exp(b*x));
end
t=(min(x):0.01:max(x));
plot(x,y,'.',t,exp(b*t),'-');
