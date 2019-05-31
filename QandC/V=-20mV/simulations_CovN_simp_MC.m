    
function [nNaT,nKT,MnNa,MnK,C0Na,C0K,QNa,QK,mINa,mIK,sINa,sIK] = simulations_CovN_simp_MC()

% coded by Roberto Fernández Galán, rfgalan@case.edu, December 2011

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   gammaNa  :   single Na channel conductance  (20 pS)
%
%   gammaK   :   single Na channel conductance  (20 pS)
%
%   memS     :   membrane patch size (example: 100 um^2)
%
%   v        :   holding potential in voltage-clamp (example value: -10 mV)
%
%   T        :   time window for integration (example value: 1000 ms)
%
%   trans    :   initial interval with nonstationary dynamics to be removed
%              (should not be smaller than 50 ms and must be larger than T)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters definition

dt = 0.01;
T = 2500;
trans = 500;

v = -20;
memS = 100;

N = round(T/dt);
Ntrans = round(trans/dt);
time = (0:1:N-1)*dt;

ENa = 50;
EK = -77;

gammaNa = 20e-9; % in mS per channel
gammaK = 20e-9;  % in mS per channel

% input parameters check

if trans >  T
    
    disp('T must be larger than the nonstationary transient\n')
    disp('Setting transient to 10%T\n')

    trans = round(T/trans);
    
end

% initializing variables

NNa = 250*memS; 
NK = 50*memS; 

% homogeneous occupancy at the beginning of the simulations

    temp = NNa/8*ones(8,1); 
    nNa = round(temp);
    [~,k] = max(nNa);
    nNa(k) = nNa(k) - (sum(nNa) - NNa);

    temp = NK/5*ones(5,1);
    nK = round(temp);
    [~,k] = max(nK);
    nK(k) = nK(k) - (sum(nK) - NK);

nNaT = zeros(8,T);
nKT = zeros(5,T);

nNaT(:,1) = nNa;
nKT(:,1) = nK;

% time integration

for t = 1:N-1

    ma = malfa(v);
    mb = mbeta(v);
    
    ha = halfa(v);
    hb = hbeta(v);
    
    na = nalfa(v);
    nb = nbeta(v);
    
    % stochastic gating of sodium channels
    
    check = 0;
    
    while check == 0
          
      PNaij = [0     nNa(1)*3*ma*dt    0       0      nNa(1)*ha*dt      0       0      0;...
               nNa(2)*mb*dt    0       nNa(2)*2*ma*dt    0      0       nNa(2)*ha*dt      0      0;...
               0     nNa(3)*2*mb*dt    0       nNa(3)*ma*dt     0       0       nNa(3)*ha*dt     0;...
               0     0       nNa(4)*3*mb*dt    0      0       0       0      sum(rand(1,round(nNa(4))) < ha*dt);...
               nNa(5)*hb*dt    0       0       0      0       nNa(5)*3*ma*dt    0      0;...
               0     nNa(6)*hb*dt      0       0      nNa(6)*mb*dt      0       nNa(6)*2*ma*dt   0;...
               0     0       nNa(7)*hb*dt      0      0       nNa(7)*2*mb*dt    0      sum(rand(1,round(nNa(7))) < ma*dt);...
               0     0       0       sum(rand(1,nNa(8)) < hb*dt)     0       0       sum(rand(1,nNa(8)) < 3*mb*dt)   0];
       
      dnNa = (PNaij' - PNaij)*ones(8,1);

      if sum(nNa + dnNa < 0) > 0
       
           check = 0;

       else

           check = 1;
           nNa = nNa + dnNa;

      end
       
      nNaT(:,t+1) = nNa;
      
    end
    
   % stochastic gating of potassium channels
    
   check = 0;
    
   while check == 0
          
      PKij = [0            nK(1)*4*na*dt   0               0                    0;...
              nK(2)*nb*dt  0               nK(2)*3*na*dt   0                    0;...
              0            nK(3)*2*nb*dt   0               nK(3)*2*na*dt        0;...
              0            0               nK(4)*3*nb*dt   0                    sum(rand(1,round(nK(4))) < na*dt);...
              0            0               0               sum(rand(1,nK(5)) < 4*nb*dt)   0];
 
       
      dnK = (PKij'-PKij)*ones(5,1);

      if sum(nK + dnK < 0) > 0
       
           check = 0;

       else

           check = 1;
           nK = nK + dnK;

      end
       
      nKT(:,t+1) = nK; 
      
   end
             
end

%% plotting traces of the occupation numbers

time = time - time(Ntrans+1);

figure
plot(time(Ntrans+1:end),nNaT(:,Ntrans+1:end))

figure
plot(time(Ntrans+1:end),nKT(:,Ntrans+1:end))

%% mean and covariances of the occupation numbers and their change

MnNa = mean(nNaT(:,Ntrans+1:N),2);
MnK = mean(nKT(:,Ntrans+1:N),2);

C0Na = cov(nNaT(:,Ntrans+1:N)');
C0K = cov(nKT(:,Ntrans+1:N)');

QNa = cov(diff(nNaT(:,Ntrans+1:N),1,2)');
QK = cov(diff(nKT(:,Ntrans+1:N),1,2)');

%% mean and standard deviation of the currents (in uA)

mINa = gammaNa*MnNa(8)*(v-ENa);
mIK = gammaK*MnK(5)*(v-EK);

sINa = gammaNa*sqrt(C0Na(8,8))*(v-ENa);
sIK = gammaK*sqrt(C0K(5,5))*(v-EK);

%% saving results

save simulations_CovN_simp_MC.mat

%% plotting means and covariances

figure
subplot(1,2,1)
bar(MnNa)
title('mean state occupancy for Na')
subplot(1,2,2)
bar(MnK)
title('mean state occupancy for K')

figure
subplot(2,2,1)
imagesc(QNa)
colorbar
title('Q_{Na}')

subplot(2,2,2)
imagesc(QK)
colorbar
title('Q_K')

subplot(2,2,3)
imagesc(C0Na)
colorbar
title('C0_{Na}')

subplot(2,2,4)
imagesc(C0K)
colorbar
title('CO_{K}')

end

%% functions for gating variables (transition rates)

function ma = malfa(V)
    
    if abs(V+40) < 0.1
		ma = 1;
	else
		ma = 0.1*(V+40)/(1-exp(-(V+40)/10));
    end
    
end

function mb = mbeta(V)

	mb = 4*exp(-(V+65)/18);

end
    
function ha = halfa(V)
	
    ha = 0.07*exp(-(V+65)/20);

end

function hb = hbeta(V)

	hb = 1/(1+exp(-(V+35)/10));

end    

function na = nalfa(V)

    if abs(V+55) < 0.1

		na = 0.1;

	else

		na = 0.01*(V+55)/(1-exp(-0.1*(V+55)));

    end
    
end
    
function nb = nbeta(V)

	nb = 0.125*exp(-(V+65)/80);
    
end
