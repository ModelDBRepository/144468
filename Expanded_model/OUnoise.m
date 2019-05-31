function noise=OUnoise(T, dt, tau, sigma)

%this function calculates Ornstein-Uhlenbeack Noise

%tau is in milliseconds

%output is in picoamps

tic

N=T/dt;
noise=zeros(1, N);
noise(1)=randn(1);

for n=2:N
    noise(n)=noise(n-1)*exp(-dt/tau)+sqrt((tau/2)*(1-exp(-2*dt/tau)))*randn(1);
end

noise=noise-mean(noise);

noise=noise/std(noise); %makes standard deviation 1

noise=noise.*sigma;

toc