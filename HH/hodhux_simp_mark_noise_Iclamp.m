function [v, I_total, I_na, I_k, timetrack, Na, K, NNa, NK,computing_time]= hodhux_simp_mark_noise_Iclamp(dt,T,Nanoise,Knoise,Nanum,Knum,trans,area)

% coded by Nicolaus Schmandt, n.schmandt@gmail.com,
% and Roberto F. Galán, rfgalan@case.edu, April 2011

%Hodgkin Huxley neuron using simplified Markov chain stochasticity


%OUTPUT


%v is the membrane trace, in millivolts

%I_total is the total current flowing through the membrane at each point in
%time, in microamps

%I_na is the total sodium current flowing at each point in time, in
%microamps

%I_k is the total potassium current flowing at each point in time, in
%microamps

%timetrack is a time index, representing time that has passed in seconds

%Na is the number of channels in each of the possible markov states

%K is the number of channels in each of the possible markov states

%NNa is the total number of sodium channels

%NK is the total number of potassium channels

%computing_time is the time the simulation took to run, in seconds


%INPUT


%dt is the timestep for integration, in milliseconds, usually between .001
%and .01

%T is the total length of the trace, in milliseconds

%Nanoise is 1 for stochastic sodium dynamics, 0 for deterministic sodium
%dynamics

%Knoise is 1 for stochastic potassium dynamics, 0 for deterministic
%potassium dynamics

%Nanum is the channel density of sodium channels, in channels per micron
%squared

%Knum is the channel density of potassium channels, in channels per micron
%squared

%trans represents "transient", an optional variable that represents how 
%long the simulation will run a deterministic transient before switching to 
%stochastic. Note that the deterministic component will still be in the 
%output.  trans is in milliseconds, if you don't want a transient, just put
%zero.

%area is the size of the neuron or membrane patch, in microns squared

tic %timing variable

%These values are the capacitance and leak constants

C = 1*(10^-8)*area; %microF/area
El = -54.387; %mVs
gl = .3*(10^-8)*area; %mS/area

%reversal potential for sodium and potassium

Ek = -77; %mVs
Ena = 45; %mVs

%sodium and potassium per channel conductance

Na_chang=20*(10^-9); %mS/channel
K_chang=20*(10^-9); %mS/channel

%total number of sodium and potassium channels

NNa=round(Nanum*area); %per area
NK=round(Knum*area); %per area

%Number of integration steps and timetrack calculation

N = round(T/dt); %time is in msec
timetrack=(1:N)*dt./1000; %timetrack is in seconds

%initializing conditions

v = zeros(N,1);

I_na=zeros(1, N);
I_k=zeros(1, N);
I_l=zeros(1, N);

I_total=zeros(1, N);

%voltage for the simulation beginning

v(1)=-60;

%calculate the starting numbers of channels in each state

K=zeros(5, N);
K(5, 1)=round(NK*(nalfa(v(1))/(nalfa(v(1))+nbeta(v(1))))^4);
K(1, 1)=NK-K(5, 1);

Na=zeros(8, N);
Na(8, 1)=round(NNa*((malfa(v(1))/(malfa(v(1))+mbeta(v(1))))^3)*(halfa(v(1))/(halfa(v(1))+hbeta(v(1)))));
Na(1, 1)=NNa-Na(8, 1);

for t = 1:N-1
    
    %integration timestep
    
    v(t+1) = v(t) + dt*(1/C)*(-I_total(t));

    %calculate alphas and betas
    
    ma=malfa(v(t));
    mb=mbeta(v(t));
    ha=halfa(v(t));
    hb=hbeta(v(t));
    na=nalfa(v(t));
    nb=nbeta(v(t));
    
    if Nanoise==0 || (t*dt)<trans %deterministic in the transient or if 
        %Nanoise is zero.
        
        %Determine transition numbers
        
        Na00_Na01=(Na(1, t)*ha*dt);
        Na01_Na00=(Na(5, t)*hb*dt);

        Na10_Na11=(Na(2, t)*ha*dt);
        Na11_Na10=(Na(6, t)*hb*dt);

        Na20_Na21=(Na(3, t)*ha*dt);
        Na21_Na20=(Na(7, t)*hb*dt);
        Na30_Na31=(Na(4, t)*ha*dt);
        Na31_Na30=(Na(8, t)*hb*dt);
        
        %Update population numbers

        Na(1, t+1)=Na(1, t)+Na01_Na00-Na00_Na01;
        Na(2, t+1)=Na(2, t)+Na11_Na10-Na10_Na11;
        Na(3, t+1)=Na(3, t)+Na21_Na20-Na20_Na21;
        Na(4, t+1)=Na(4, t)+Na31_Na30-Na30_Na31;

        Na(5, t+1)=Na(5, t)+Na00_Na01-Na01_Na00;
        Na(6, t+1)=Na(6, t)+Na10_Na11-Na11_Na10;
        Na(7, t+1)=Na(7, t)+Na20_Na21-Na21_Na20;
        Na(8, t+1)=Na(8, t)+Na30_Na31-Na31_Na30;
        
        %Determine transition numbers

        Na00_Na10=(Na(1, t+1)*3*ma*dt);
        Na10_Na00=(Na(2, t+1)*mb*dt);

        Na10_Na20=(Na(2, t+1)*2*ma*dt);
        Na20_Na10=(Na(3, t+1)*2*mb*dt);

        Na20_Na30=(Na(3, t+1)*ma*dt);
        Na30_Na20=(Na(4, t+1)*3*mb*dt);

        Na01_Na11=(Na(5, t+1)*3*ma*dt);
        Na11_Na01=(Na(6, t+1)*mb*dt);

        Na11_Na21=(Na(6, t+1)*2*ma*dt);
        Na21_Na11=(Na(7, t+1)*2*mb*dt);
        Na21_Na31=(Na(7, t+1)*ma*dt);
        Na31_Na21=(Na(8, t+1)*3*mb*dt);
        
        %Update population numbers

        Na(1, t+1)=Na(1, t+1)+Na10_Na00-Na00_Na10;
        Na(2, t+1)=Na(2, t+1)+Na00_Na10+Na20_Na10-Na10_Na20-Na10_Na00;
        Na(3, t+1)=Na(3, t+1)+Na10_Na20+Na30_Na20-Na20_Na30-Na20_Na10;
        Na(4, t+1)=Na(4, t+1)+Na20_Na30-Na30_Na20;

        Na(5, t+1)=Na(5, t+1)+Na11_Na01-Na01_Na11;
        Na(6, t+1)=Na(6, t+1)+Na01_Na11+Na21_Na11-Na11_Na21-Na11_Na01;
        Na(7, t+1)=Na(7, t+1)+Na11_Na21+Na31_Na21-Na21_Na11-Na21_Na31;
        Na(8, t+1)=Na(8, t+1)+Na21_Na31-Na31_Na21;

        temp=1; 

        while Na(1, t+1)<0 || Na(2, t+1)<0 || Na(3, t+1)<0 || Na(4, t+1)<0 ...
                || Na(5, t+1)<0 || Na(6, t+1)<0 || Na(7, t+1)<0 || Na(8, t+1)<0 %check to make sure the number of channels in all states is positive
            
            %Determine transition numbers

            Na00_Na01=(Na(1, t)*ha*dt)/temp;
            Na01_Na00=(Na(5, t)*hb*dt)/temp;

            Na10_Na11=(Na(2, t)*ha*dt)/temp;
            Na11_Na10=(Na(6, t)*hb*dt)/temp;

            Na20_Na21=(Na(3, t)*ha*dt)/temp;
            Na21_Na20=(Na(7, t)*hb*dt)/temp;

            Na30_Na31=(Na(4, t)*ha*dt)/temp;
            Na31_Na30=(Na(8, t)*hb*dt)/temp;

            %Update population numbers

            Na(1, t+1)=Na(1, t)+Na01_Na00-Na00_Na01;
            Na(2, t+1)=Na(2, t)+Na11_Na10-Na10_Na11;
            Na(3, t+1)=Na(3, t)+Na21_Na20-Na20_Na21;
            Na(4, t+1)=Na(4, t)+Na31_Na30-Na30_Na31;

            Na(5, t+1)=Na(5, t)+Na00_Na01-Na01_Na00;
            Na(6, t+1)=Na(6, t)+Na10_Na11-Na11_Na10;
            Na(7, t+1)=Na(7, t)+Na20_Na21-Na21_Na20;
            Na(8, t+1)=Na(8, t)+Na30_Na31-Na31_Na30;
            
            %Determine transition numbers

            Na00_Na10=(Na(1, t+1)*3*ma*dt)/temp;
            Na10_Na00=(Na(2, t+1)*mb*dt)/temp;

            Na10_Na20=(Na(2, t+1)*2*ma*dt)/temp;
            Na20_Na10=(Na(3, t+1)*2*mb*dt)/temp;

            Na20_Na30=(Na(3, t+1)*ma*dt)/temp;
            Na30_Na20=(Na(4, t+1)*3*mb*dt)/temp;


            Na01_Na11=(Na(5, t+1)*3*ma*dt)/temp;
            Na11_Na01=(Na(6, t+1)*mb*dt)/temp;

            Na11_Na21=(Na(6, t+1)*2*ma*dt)/temp;
            Na21_Na11=(Na(7, t+1)*2*mb*dt)/temp;

            Na21_Na31=(Na(7, t+1)*ma*dt)/temp;
            Na31_Na21=(Na(8, t+1)*3*mb*dt)/temp;

            %Update population numbers

            Na(1, t+1)=Na(1, t+1)+Na10_Na00-Na00_Na10;
            Na(2, t+1)=Na(2, t+1)+Na00_Na10+Na20_Na10-Na10_Na20-Na10_Na00;
            Na(3, t+1)=Na(3, t+1)+Na10_Na20+Na30_Na20-Na20_Na30-Na20_Na10;
            Na(4, t+1)=Na(4, t+1)+Na20_Na30-Na30_Na20;

            Na(5, t+1)=Na(5, t+1)+Na11_Na01-Na01_Na11;
            Na(6, t+1)=Na(6, t+1)+Na01_Na11+Na21_Na11-Na11_Na21-Na11_Na01;
            Na(7, t+1)=Na(7, t+1)+Na11_Na21+Na31_Na21-Na21_Na11-Na21_Na31;
            Na(8, t+1)=Na(8, t+1)+Na21_Na31-Na31_Na21;

            temp=temp+1; %this reduces the amount of each change in the event of one of the channel numbers going negative.  Everytime I have checked, this loop has never been used, but I leave it in for consistency.

        end
        
    else
        
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
        
    end
    
    if Knoise==0 || (t*dt)<trans %deterministic in the transient or if 
        %Knoise is zero.
        
        %Determine transition numbers
        
        K1_K2=(K(1, t)*4*na*dt);
        K2_K1=(K(2, t)*nb*dt);

        K2_K3=(K(2, t)*3*na*dt);
        K3_K2=(K(3, t)*2*nb*dt);

        K3_K4=(K(3, t)*2*na*dt);
        K4_K3=(K(4, t)*3*nb*dt);

        K4_K5=K(4, t)*na*dt;
        K5_K4=K(5, t)*4*nb*dt;
        
        %Update channel populations

        K(1, t+1)=K(1, t)+K2_K1-K1_K2;
        K(2, t+1)=K(2, t)+K1_K2+K3_K2-K2_K1-K2_K3;
        K(3, t+1)=K(3, t)+K2_K3+K4_K3-K3_K2-K3_K4;
        K(4, t+1)=K(4, t)+K3_K4+K5_K4-K4_K3-K4_K5;
        K(5, t+1)=K(5, t)+K4_K5-K5_K4;

        temp=1;

        while K(1, t+1)<0 || K(2, t+1)<0 || K(3, t+1)<0 || K(4, t+1)<0 || K(5, t+1)<0 %check to make sure the number of channels in all states is positive
            
            %Determine transition numbers

            K1_K2=(K(1, t)*4*na*dt)/temp;
            K2_K1=(K(2, t)*nb*dt)/temp;

            K2_K3=(K(2, t)*3*na*dt)/temp;
            K3_K2=(K(3, t)*2*nb*dt)/temp;

            K3_K4=(K(3, t)*2*na*dt)/temp;
            K4_K3=(K(4, t)*3*nb*dt)/temp;

            K4_K5=K(4, t)*na*dt/temp;
            K5_K4=K(5, t)*4*nb*dt/temp;
            
            %Update channel populations

            K(1, t+1)=K(1, t)+K2_K1-K1_K2;
            K(2, t+1)=K(2, t)+K1_K2+K3_K2-K2_K1-K2_K3;
            K(3, t+1)=K(3, t)+K2_K3+K4_K3-K3_K2-K3_K4;
            K(4, t+1)=K(4, t)+K3_K4+K5_K4-K4_K3-K4_K5;
            K(5, t+1)=K(5, t)+K4_K5-K5_K4;

            temp=temp+1;

        end
        
    else
        
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
        
    end
    
    %calculate currents at each time step
    
    I_na(t+1)=Na(8, t+1)*Na_chang*(v(t+1)-Ena);
    I_k(t+1)=K(5, t+1)*K_chang*(v(t+1)-Ek);
    I_l(t+1)=gl*(v(t+1)-El);
    
    %calculate total current at each time step
    
    I_total(t+1)=I_na(t+1)+I_k(t+1)+I_l(t+1);
        
end    

computing_time=toc; %timing variable
toc
        
function ma = malfa(V) %From Dayan and Abbot textbook
    
    if abs(V+40)<0.1
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

    if abs(V+55)<0.1

		na = 0.1;

	else

		na = 0.01*(V+55)/(1-exp(-0.1*(V+55)));

    end
    
end
    

function nb = nbeta(V) %From Dayan and Abbot textbook

	nb = 0.125*exp(-(V+65)/80);
    
end

end