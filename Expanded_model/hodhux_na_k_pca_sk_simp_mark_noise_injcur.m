function [v, I_na, I_k, I_pca, I_sk, I_l, I_total, caconc, timetrack, NP, Na, K, PCa, SK, g_current]= hodhux_na_k_pca_sk_simp_mark_noise_injcur(dt,T,Nanum,Knum,Pnum,SKnum,Idc,g_input,buffering_constant,catau,area)

%Coded by Nicolaus Schmandt, n.schmandt@gmail.com, May 2012
%
%---input variables---
%
%dt=timestep in milliseconds (usually somewhere between .001 and .01
%
%T=total run time in milliseconds
%
%Nanoise,Knoise,PCanoise,SKnoise=a binary variable that
%determines whether each current is stochastic or deterministic, with 1
%being stochastic and 0 deterministic
%
%Nanum,Knum,Pnum,SKnum=channel densities of each of the channels
%for the different kinds of currents.
%
%trans=a deterministic transient where all currents will be deterministic
%before switching to stochastic. This will be included in the output
%
%Idc=injected or received currents, as a vector of length equal to the
%number of time steps.
%
%g_input=input synaptic conductances (assumed to have a reversal potential
%of zero) in mS
%
%buffering_constant=representation of the amount of calcium buffering in a
%neuron, the higher the buffering, the less calcium per spike will affect
%the total calcium concentration.
%
%catau=speed with which calcium is removed from the cell, as a single
%exponential decay.
%
%area=size of the neuron
%
%---output variables---
%
%v=vector of membrane potential output
%
%I_na,I_k,I_pca,I_sk,I_l,I_total=sodium, potassium, p-type calcium, 
%sk-potassium, leak, and total current vectors
%
%caconc=calcium concentration in the neuron
%
%timetrack=vector of time steps, in seconds
%
%NP=number of P-type calcium channels
%
%Na, K, PCa, SK=number of each channels in each of the markov
%chain states
%
%g_current=current due to synaptic input

tic

% parameter definition

C = 1*(10^-8)*area; %microF/area
El = -54.387; %mVs
gl = .3*(10^-8)*area; %mS/area

Ek = -77; %mVs
Ena = 45; %mVs

Na_chang=20*(10^-9); %mS/channel
K_chang=20*(10^-9); %mS/channel
P_perm=.3; %microns/s based on Stanley et. al.'s value for Cav2.2 (N-type) 2010, per channel per micron^2
SK_chang=10*(10^-9); %mS/channel, from Faber and Sah 2007

NNa=Nanum*area; %per area
NK=Knum*area; %per area
NP=round(Pnum*area); %per area
NSK=round(SKnum*area); 

N = T/dt; %time is in msec
timetrack=(1:N)*dt./1000;

v=zeros(N,1);

basal_calcium=.1;

caconc=zeros(1, N);
caconc(1) = basal_calcium; %caconc is in micromoles

v(1)=-60;

pca_starting_fraction=.9;

Na=zeros(8, N);
Na(8, 1)=round(NNa*((malfa(v(1))/(malfa(v(1))+mbeta(v(1))))^3)*(halfa(v(1))/(halfa(v(1))+hbeta(v(1)))));
Na(1, 1)=NNa-Na(8, 1);

K=zeros(5, N);
K(5, 1)=round(NK*((nalfa(v(1))/(nalfa(v(1))+nbeta(v(1))))^4));
K(1, 1)=NK-K(5, 1);

PCa=zeros(4, N);
PCa(3, 1)=round(pca_starting_fraction*NP);
PCa(1, 1)=NP-PCa(3, 1);

SK=zeros(6, N);

SK(1, 1)=NSK;

% initializing the dynamics

I_na=zeros(1, N);
I_k=zeros(1, N);
I_l=zeros(1, N);
I_pca=zeros(1, N);
I_sk=zeros(1, N);

ca_current=zeros(1, N);
k_current=zeros(1, N);

p_ca_current=zeros(1, N);

p_total_current=zeros(1, N);

I_total=zeros(1, N);

%assuming radius of 3 microns
volume=buffering_constant*(4/3)*pi*(sqrt(area/(4*pi)))^3; %microns^3

gamma=(1/(2*96485))*(1/1000)*(10^15)/volume; %1/Faraday's Constant, a correction to milliseconds, correction from microns^3 to Liters, divided by the volume

ph_tau=1500; %in milliseconds

p_ca_ratio=1/5000; %ratio of calcium to potassium conductance in p-channels

ca_current(1)=(1/(10^15))*(4*v(1)*96485/26.7)*(caconc(1)-1300*exp(-2*v(1)/26.7))/(1-exp(-2*v(1)/26.7)); %in microamps, based on GHK current theory
%(conversion from Liters to microns^3, z^2, voltage, Faraday's
%constant/(RT/F))*(internal caconc(in micromolar)-external caconc*exp)
k_current(1)=(1/(10^15))*(4*v(1)*96485/26.7)*(155000-4700*exp(-2*v(1)/26.7))/(1-exp(-2*v(1)/26.7)); %in microamps

if v(1)==0
    ca_current(1)=2*(96485)*(caconc(1)-1300)/(10^15);
    k_current(1)=2*(96485)*(155000-4700)/(10^15);
end

p_ca_current(1)=P_perm*ca_current(1);
p_total_current(1)=p_ca_current(1)+p_ca_ratio*k_current(1);

g_current=zeros(1, N);

syn_reversal=0; %mV

for t = 1:N-1
    
    v(t+1) = v(t) + dt*(1/C)*(-I_total(t)+Idc(t)-g_input(t)*(v(t)-syn_reversal));
    
    g_current(t+1)=g_input(t)*(v(t)-syn_reversal);
    
    ma=malfa(v(t));
    mb=mbeta(v(t));
    ha=halfa(v(t));
    hb=hbeta(v(t));
    na=nalfa(v(t));
    nb=nbeta(v(t));
    
    pm_inf=pm_infcalc(v(t));
    ph_inf=ph_infcalc(v(t));
    pm_tau=pm_taucalc(v(t));
    
    pcama=pm_inf/pm_tau;
    pcamb=(1-pm_inf)/pm_tau;
    pcaha=ph_inf/ph_tau;
    pcahb=(1-ph_inf)/ph_tau;
    
    ca_current(t+1)=(1/(10^15))*(4*v(t)*96485/26.7)*(caconc(t)-1300*exp(-2*v(t)/26.7))/(1-exp(-2*v(t)/26.7)); %in microamps, based on GHK current theory
    %(conversion from Liters to microns^3, z^2, voltage, Faraday's
    %constant/(RT/F))*(internal caconc(in micromolar)-external caconc*exp)
    k_current(t+1)=(1/(10^15))*(4*v(t)*96485/26.7)*(155000-4700*exp(-2*v(t)/26.7))/(1-exp(-2*v(t)/26.7)); %in microamps
    
    if v(t)==0
        ca_current(t+1)=2*(96485)*(caconc(t)-1300)/(10^15);
        k_current(t+1)=2*(96485)*(155000-4700)/(10^15);
    end
    
    p_ca_current(t+1)=P_perm*ca_current(t+1);
    p_total_current(t+1)=p_ca_current(t+1)+p_ca_ratio*k_current(t+1);
    
    Na1_Na5=(Na(1, t)*ha*dt);
    Na5_Na1=(Na(5, t)*hb*dt);

    Na2_Na6=(Na(2, t)*ha*dt);
    Na6_Na2=(Na(6, t)*hb*dt);

    Na3_Na7=(Na(3, t)*ha*dt);
    Na7_Na3=(Na(7, t)*hb*dt);

    Na4_Na8=sum(rand(round(Na(4, t)), 1)<ha*dt);
    Na8_Na4=sum(rand(round(Na(8, t)), 1)<hb*dt);

    Na1_Na2=(Na(1, t)*3*ma*dt);
    Na2_Na1=(Na(2, t)*mb*dt);

    Na2_Na3=(Na(2, t)*2*ma*dt);
    Na3_Na2=(Na(3, t)*2*mb*dt);

    Na3_Na4=(Na(3, t)*ma*dt);
    Na4_Na3=(Na(4, t)*3*mb*dt);

    Na5_Na6=(Na(5, t)*3*ma*dt);
    Na6_Na5=(Na(6, t)*mb*dt);

    Na6_Na7=(Na(6, t)*2*ma*dt);
    Na7_Na6=(Na(7, t)*2*mb*dt);

    Na7_Na8=sum(rand(round(Na(7, t)), 1)<ma*dt);
    Na8_Na7=sum(rand(round(Na(8, t)), 1)<3*mb*dt);

    %Update channel populations

    Na(1, t+1)=Na(1, t)+Na5_Na1+Na2_Na1-Na1_Na2-Na1_Na5;
    Na(2, t+1)=Na(2, t)+Na6_Na2+Na1_Na2+Na3_Na2-Na2_Na3-Na2_Na1-Na2_Na6;
    Na(3, t+1)=Na(3, t)+Na7_Na3+Na2_Na3+Na4_Na3-Na3_Na4-Na3_Na2-Na3_Na7;
    Na(4, t+1)=Na(4, t)+Na8_Na4+Na3_Na4-Na4_Na3-Na4_Na8;

    Na(5, t+1)=Na(5, t)+Na1_Na5+Na6_Na5-Na5_Na6-Na5_Na1;
    Na(6, t+1)=Na(6, t)+Na2_Na6+Na5_Na6+Na7_Na6-Na6_Na7-Na6_Na5-Na6_Na2;
    Na(7, t+1)=Na(7, t)+Na3_Na7+Na6_Na7+Na8_Na7-Na7_Na6-Na7_Na8-Na7_Na3;
    Na(8, t+1)=Na(8, t)+Na4_Na8+Na7_Na8-Na8_Na7-Na8_Na4;

    while Na(1, t+1)<0 || Na(2, t+1)<0 || Na(3, t+1)<0 || Na(4, t+1)<0 ...
            || Na(5, t+1)<0 || Na(6, t+1)<0 || Na(7, t+1)<0 || Na(8, t+1)<0 %check to make sure the number of channels in all states is positive

        %Determine transition numbers

        Na1_Na5=(Na(1, t)*ha*dt);
        Na5_Na1=(Na(5, t)*hb*dt);

        Na2_Na6=(Na(2, t)*ha*dt);
        Na6_Na2=(Na(6, t)*hb*dt);

        Na3_Na7=(Na(3, t)*ha*dt);
        Na7_Na3=(Na(7, t)*hb*dt);

        Na4_Na8=sum(rand(round(Na(4, t)), 1)<ha*dt);
        Na8_Na4=sum(rand(round(Na(8, t)), 1)<hb*dt);

        Na1_Na2=(Na(1, t)*3*ma*dt);
        Na2_Na1=(Na(2, t)*mb*dt);

        Na2_Na3=(Na(2, t)*2*ma*dt);
        Na3_Na2=(Na(3, t)*2*mb*dt);

        Na3_Na4=(Na(3, t)*ma*dt);
        Na4_Na3=(Na(4, t)*3*mb*dt);

        Na5_Na6=(Na(5, t)*3*ma*dt);
        Na6_Na5=(Na(6, t)*mb*dt);

        Na6_Na7=(Na(6, t)*2*ma*dt);
        Na7_Na6=(Na(7, t)*2*mb*dt);

        Na7_Na8=sum(rand(round(Na(7, t)), 1)<ma*dt);
        Na8_Na7=sum(rand(round(Na(8, t)), 1)<3*mb*dt);

        %Update channel populations

        Na(1, t+1)=Na(1, t)+Na5_Na1+Na2_Na1-Na1_Na2-Na1_Na5;
        Na(2, t+1)=Na(2, t)+Na6_Na2+Na1_Na2+Na3_Na2-Na2_Na3-Na2_Na1-Na2_Na6;
        Na(3, t+1)=Na(3, t)+Na7_Na3+Na2_Na3+Na4_Na3-Na3_Na4-Na3_Na2-Na3_Na7;
        Na(4, t+1)=Na(4, t)+Na8_Na4+Na3_Na4-Na4_Na3-Na4_Na8;

        Na(5, t+1)=Na(5, t)+Na1_Na5+Na6_Na5-Na5_Na6-Na5_Na1;
        Na(6, t+1)=Na(6, t)+Na2_Na6+Na5_Na6+Na7_Na6-Na6_Na7-Na6_Na5-Na6_Na2;
        Na(7, t+1)=Na(7, t)+Na3_Na7+Na6_Na7+Na8_Na7-Na7_Na6-Na7_Na8-Na7_Na3;
        Na(8, t+1)=Na(8, t)+Na4_Na8+Na7_Na8-Na8_Na7-Na8_Na4;

    end

    K1_K2=(K(1, t)*4*na*dt);
    K2_K1=(K(2, t)*nb*dt);

    K2_K3=(K(2, t)*3*na*dt);
    K3_K2=(K(3, t)*2*nb*dt);

    K3_K4=(K(3, t)*2*na*dt);
    K4_K3=(K(4, t)*3*nb*dt);

    K4_K5=sum(rand(round(K(4, t)), 1)<na*dt);
    K5_K4=sum(rand(round(K(5, t)), 1)<4*nb*dt);

    %Update channel populations

    K(1, t+1)=K(1, t)+K2_K1-K1_K2;
    K(2, t+1)=K(2, t)+K1_K2+K3_K2-K2_K1-K2_K3;
    K(3, t+1)=K(3, t)+K2_K3+K4_K3-K3_K2-K3_K4;
    K(4, t+1)=K(4, t)+K3_K4+K5_K4-K4_K3-K4_K5;
    K(5, t+1)=K(5, t)+K4_K5-K5_K4;

    while K(1, t+1)<0 || K(2, t+1)<0 || K(3, t+1)<0 || K(4, t+1)<0 || K(5, t+1)<0 %check to make sure the number of channels in all states is positive

        %Determine transition numbers

        K1_K2=(K(1, t)*4*na*dt);
        K2_K1=(K(2, t)*nb*dt);

        K2_K3=(K(2, t)*3*na*dt);
        K3_K2=(K(3, t)*2*nb*dt);

        K3_K4=(K(3, t)*2*na*dt);
        K4_K3=(K(4, t)*3*nb*dt);

        K4_K5=sum(rand(round(K(4, t)), 1)<na*dt);
        K5_K4=sum(rand(round(K(5, t)), 1)<4*nb*dt);

        %Update channel populations

        K(1, t+1)=K(1, t)+K2_K1-K1_K2;
        K(2, t+1)=K(2, t)+K1_K2+K3_K2-K2_K1-K2_K3;
        K(3, t+1)=K(3, t)+K2_K3+K4_K3-K3_K2-K3_K4;
        K(4, t+1)=K(4, t)+K3_K4+K5_K4-K4_K3-K4_K5;
        K(5, t+1)=K(5, t)+K4_K5-K5_K4;

    end
    
    P1_P2=PCa(1,t)*pcama*dt;
    P2_P1=PCa(2,t)*pcamb*dt;
    
    P1_P3=PCa(1,t)*pcaha*dt;
    P3_P1=PCa(3,t)*pcahb*dt;
    
    P2_P4=sum(rand(round(PCa(2,t)),1)<pcaha*dt);
    P4_P2=sum(rand(round(PCa(4,t)),1)<pcahb*dt);
    
    P3_P4=sum(rand(round(PCa(3,t)),1)<pcama*dt);
    P4_P3=sum(rand(round(PCa(4,t)),1)<pcamb*dt);

    PCa(1, t+1)=PCa(1, t)+P2_P1+P3_P1-P1_P2-P1_P3;
    PCa(2, t+1)=PCa(2, t)+P1_P2+P4_P2-P2_P4-P2_P1;
    PCa(3, t+1)=PCa(3, t)+P1_P3+P4_P3-P3_P1-P3_P4;
    PCa(4, t+1)=PCa(4, t)+P2_P4+P3_P4-P4_P2-P4_P3;
        
    SK1_SK2=SK(1, t)*.2*dt*caconc(t);
    SK2_SK1=SK(2, t)*.08*dt;

    SK2_SK3=SK(2, t)*.16*dt*caconc(t);
    SK3_SK2=SK(3, t)*.08*dt;

    SK3_SK4=SK(3, t)*.08*dt*caconc(t);
    SK4_SK3=SK(4, t)*.2*dt;

    SK3_SK5=sum(rand(round(SK(3, t)), 1)<.16*dt);
    SK5_SK3=sum(rand(round(SK(5, t)), 1)<1*dt);

    SK4_SK6=sum(rand(round(SK(4, t)), 1)<1.2*dt);
    SK6_SK4=sum(rand(round(SK(6, t)), 1)<.1*dt);

    SK(1, t+1)=SK(1, t)+SK2_SK1-SK1_SK2;
    SK(2, t+1)=SK(2, t)+SK1_SK2+SK3_SK2-SK2_SK1-SK2_SK3;
    SK(3, t+1)=SK(3, t)+SK2_SK3+SK5_SK3+SK4_SK3-SK3_SK2-SK3_SK5-SK3_SK4;
    SK(4, t+1)=SK(4, t)+SK3_SK4+SK6_SK4-SK4_SK3-SK4_SK6;
    SK(5, t+1)=SK(5, t)+SK3_SK5-SK5_SK3;
    SK(6, t+1)=SK(6, t)+SK4_SK6-SK6_SK4;

    while SK(1, t+1)<0 || SK(2, t+1)<0 || SK(3, t+1)<0 || SK(4, t+1)<0 || SK(5, t+1)<0 || SK(6, t+1)<0

        SK1_SK2=SK(1, t)*.2*dt*caconc(t);
        SK2_SK1=SK(2, t)*.08*dt;

        SK2_SK3=SK(2, t)*.16*dt*caconc(t);
        SK3_SK2=SK(3, t)*.08*dt;

        SK3_SK4=SK(3, t)*.08*dt*caconc(t);
        SK4_SK3=SK(4, t)*.2*dt;

        SK3_SK5=sum(rand(round(SK(3, t)), 1)<.16*dt);
        SK5_SK3=sum(rand(round(SK(5, t)), 1)<1*dt);

        SK4_SK6=sum(rand(round(SK(4, t)), 1)<1.2*dt);
        SK6_SK4=sum(rand(round(SK(6, t)), 1)<.1*dt);

        SK(1, t+1)=SK(1, t)+SK2_SK1-SK1_SK2;
        SK(2, t+1)=SK(2, t)+SK1_SK2+SK3_SK2-SK2_SK1-SK2_SK3;
        SK(3, t+1)=SK(3, t)+SK2_SK3+SK5_SK3+SK4_SK3-SK3_SK2-SK3_SK5-SK3_SK4;
        SK(4, t+1)=SK(4, t)+SK3_SK4+SK6_SK4-SK4_SK3-SK4_SK6;
        SK(5, t+1)=SK(5, t)+SK3_SK5-SK5_SK3;
        SK(6, t+1)=SK(6, t)+SK4_SK6-SK6_SK4;

    end
        
    dca = dt*(-(PCa(4, t)*p_ca_current(t))*gamma-((caconc(t)-basal_calcium)/catau));   
    caconc(t+1)=caconc(t)+dca; %caconc is in micromoles per liter, only in the dendrites
    
    I_na(t+1)=Na(8, t+1)*Na_chang*(v(t+1)-Ena);
    I_k(t+1)=K(5, t+1)*K_chang*(v(t+1)-Ek);
    I_l(t+1)=gl*(v(t+1)-El);
    I_pca(t+1)=PCa(4, t+1)*p_total_current(t);
    I_sk(t+1)=(SK(5, t+1)+SK(6, t+1))*SK_chang*(v(t+1)-Ek);

    I_total(t+1)=I_na(t+1)+I_k(t+1)+I_l(t+1)+I_pca(t+1)+I_sk(t+1);
        
end    

toc

%% gating variables

function ma = malfa(V) %From Dayan and Abbot textbook
    
    if abs(V+40)<1
		ma = 1;
	else
		ma = 0.1*(V+40)/(1-exp(-(V+40)/10));
    end
    
end
    

function mb = mbeta(V) %From Dayan and Abbot textbook

	mb = 4*exp(-(V+65)/18);

end
    
function ha = halfa(V) %From Dayan and Abbot textbook
	
    ha = 0.07*exp(-(V+65)/20);

end

function hb = hbeta(V) %From Dayan and Abbot textbook

	hb = 1/(1+exp(-(V+35)/10));

end    


function na = nalfa(V) %From Dayan and Abbot textbook

    if abs(V+55)<1

		na = 0.1;

	else

		na = 0.01*(V+55)/(1-exp(-0.1*(V+55)));

    end
    
end
    

function nb = nbeta(V) %From Dayan and Abbot textbook

	nb = 0.125*exp(-(V+65)/80);
    
end

function pm_inf=pm_infcalc(v) %set to match Regan Fig.7

    %pm_inf=1./(1+exp((v+9)/-5.8));
    pm_inf=1./(1+exp((v+19)/-5.8));
    
end

function ph_inf=ph_infcalc(v)

    %ph_inf=(1./(1+exp((v+21)/6.9)));
    ph_inf=(1./(1+exp((v+31)/6.9)));
    
end

function pm_tau = pm_taucalc(v)

    %pm_tau=(1/(exp((v+44.47)/-19.41)+exp((v-7.32)/13.72)))+.195; %Wakamori fig.9 and Regan fig.4
    pm_tau=(1/(exp((v+54.47)/-19.41)+exp((v+2.68)/13.72)))+.195; 
    
end

end