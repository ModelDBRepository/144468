
%This file Generates a Poisson Spike Train Output

%Specify T (milliseconds) in the command prompt to set the total time

%Specify lambda (milliseconds) in the command prompt to set the Frequency at which events occur

%The variable "cond" in the output represents the input train

%Specify the max_cond to set the peak event input value

dt=.01;
time=T;
N=time/dt;
exp_length=.5;

max_cond=10.*10^-6; 
tau=1; %msec
u=rand(N,1)<dt/lambda;
w=(dt*(1:1:10000)./tau).*exp(-dt*(1:1:10000)/tau);
u=conv(double(u),w);
u=u(1:N);
u=u/std(u);

cond=max_cond*u;