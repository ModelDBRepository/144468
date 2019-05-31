function [VRNa,VLNa,VRK,VLK,DNa,DK,QNa,QK,C0Na,C0K,ANa,AK,mINa,mIK,sINa,sIK] = theoretical_CovN_simp_MC()

% coded by Roberto Fernández Galán, rfgalan@case.edu, December 2011

clear
close all
clc

%% Parameter definition

v = -20;
dt = 0.01;
memS = 100;

ENa = 50;
EK = -77;

gammaNa = 20e-9;
gammaK = 20e-9;

NNa = 250*memS;
NK = 50*memS;

ma = malfa(v);
mb = mbeta(v);
    
ha = halfa(v);
hb = hbeta(v);
    
na = nalfa(v);
nb = nbeta(v);

%% Transition matrices

ANa = [0     3*ma    0       0      ha      0       0      0;...
       mb    0       2*ma    0      0       ha      0      0;...
       0     2*mb    0       ma     0       0       ha     0;...
       0     0       3*mb    0      0       0       0      ha;...
       hb    0       0       0      0       3*ma    0      0;...
       0     hb      0       0      mb      0       2*ma   0;...
       0     0       hb      0      0       2*mb    0      ma;...
       0     0       0       hb     0       0       3*mb   0];
   
AK = [0      4*na    0       0      0;...
      nb     0       3*na    0      0;...
      0      2*nb    0       2*na   0;...
      0      0       3*nb    0      na;...
      0      0       0       4*nb   0];  
  
  
      ANasimp = zeros(8);
      ANasimp([4 7 8],[4 7 8]) = ANa([4 7 8],[4 7 8]);

      AKsimp = zeros(5);
      AKsimp([4 5],[4 5]) = AK([4 5],[4 5]);  
  
  
ANa = ANa' - diag(sum(ANa,2));
AK  = AK'  - diag(sum(AK,2));

       ANasimp = ANasimp' - diag(sum(ANasimp,2));
       AKsimp  = AKsimp'  - diag(sum(AKsimp,2));
   
nNa = size(ANa,1);
nK = size(AK,1);
  
%% Calculation of the steady states

NNas = null(ANa);
NKs  = null(AK);

NNas = NNas/sum(NNas)*NNa;
NKs = NKs/sum(NKs)*NK;

%% Eigenvalue decomposition of the transition matrices

% Sodium

[VRNa,DNa] = eig(ANa);       
[DNa,ind] = sort(diag(DNa),'ascend');
DNa = diag(DNa);
VRNa = VRNa(:,ind);

[VLNa,DNa] = eig(ANa.');
VLNa = conj(VLNa);
[DNa,ind] = sort(diag(DNa),'ascend');
DNa = diag(DNa);
VLNa = VLNa(:,ind);

for n = 1:nNa
    
    VLNa(:,n) = VLNa(:,n)/(VLNa(:,n).'*VRNa(:,n));
    
end

% Potassium

[VRK,DK] = eig(AK);
[DK,ind] = sort(diag(DK),'ascend');
DK = diag(DK);
VRK = VRK(:,ind);

[VLK,DK] = eig(AK.');
VLK = conj(VLK);
[DK,ind] = sort(diag(DK),'ascend');
DK = diag(DK);
VLK = VLK(:,ind);

for n = 1:nK
    
    VLK(:,n) = VLK(:,n)/(VLK(:,n).'*VRK(:,n));
    
end

%% Covariance of the transition fluctuations 

% Sodium

QNa = zeros(nNa);

for i = [4 7 8]
    
    for j = [4 7 8] 
    
        if i~=j
        
            QNa(i,j) = -( NNas(j)*ANasimp(i,j) + NNas(i)*ANasimp(j,i) )*dt;   
            
        else
            
            temp1 = NNas;
            temp1(i) = [];
            
            temp2 = ANasimp(i,:);
            temp2(i) = [];
            
            QNa(i,i) = ( temp2*temp1 - NNas(i)*ANasimp(i,i) )*dt;
            
        end 
        
    end
    
end

% Potassium

QK = zeros(nK);

for i = [4 5]
    
    for j = [4 5]
    
        if i~=j
        
            QK(i,j) = -( NKs(j)*AKsimp(i,j) + NKs(i)*AKsimp(j,i) )*dt;   
            
        else
            
            temp1 = NKs;
            temp1(i) = [];
            
            temp2 = AKsimp(i,:);
            temp2(i) = [];
            
            QK(i,i) = ( temp2*temp1 - NKs(i)*AKsimp(i,i) )*dt;
            
        end  
        
    end
    
end

%% Calculating the covariance of the occupation numbers

% sodium

VRNa = VRNa(:,1:nNa-1);
VLNa = VLNa(:,1:nNa-1);

DNa = diag(DNa);
DNa = DNa(1:end-1);
C0Na = (VLNa'*QNa*VLNa)./((DNa*dt)*ones(1,nNa-1)+ones(nNa-1,1)*(DNa*dt)');
C0Na = -VRNa*C0Na*VRNa';

% potassium

VRK = VRK(:,1:nK-1);
VLK = VLK(:,1:nK-1);

DK = diag(DK);
DK = DK(1:end-1);
C0K = (VLK'*QK*VLK)./((DK*dt)*ones(1,nK-1)+ones(nK-1,1)*(DK*dt)');
C0K = -VRK*C0K*VRK';

%% mean and standard deviation of the currents (in uA)

mINa = gammaNa*NNas(8)*(v-ENa);
mIK = gammaK*NKs(5)*(v-EK);

sINa = gammaNa*sqrt(C0Na(8,8))*(v-ENa);
sIK = gammaK*sqrt(C0K(5,5))*(v-EK);

%% Saving results

save theoretical_CovN_simp_MC.mat

%% Plotting results

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
title('C0_{K}')

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