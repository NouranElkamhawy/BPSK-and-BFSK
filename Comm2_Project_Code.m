%constants
numberOfBits = 20;
Tb = 40;
periodIntervals = 1:Tb;
No = 2;
Ts = 1;
SNR = [-4:4];
Wc = 4*((2*pi)/Tb);
no = 4;
n1 = 1;
W1 = 2*pi*(no+n1)/Tb;
W2 = 2*pi*(no-n1)/Tb;
ErrorBPSK = 0;
AverageErrorBPSK = 0;
ErrorBFSK = 0;
AverageErrorBFSK = 0;
BPSK_TH_BE = 0;
BFSK_TH_BE = 0;

%-------------------------------------

% generate a random array containing 0 and 1 only 
BPSK_bits = randi ([0,1],1,numberOfBits);
BFSK_bits = BPSK_bits;
% replace each 0 in the array with one 
BPSK_bits(BPSK_bits == 0)= -1;




%generate pulses that takes time duration Tb
P_pulse = rectpulse(BPSK_bits,Tb);
F_pulse = rectpulse(BFSK_bits,Tb);

% get the effective bandwidth
[BPSK_PSD, BPSK_Freq] = periodogram( P_pulse, [],[],1);
figure
plot(BPSK_Freq,BPSK_PSD)
xlabel('Frequency') 
ylabel('Power Spectral Density')  

[BFSK_PSD, BFSK_Freq] = periodogram( F_pulse, [],[],1);
figure
plot(BFSK_Freq,BFSK_PSD)
xlabel('Frequency') 
ylabel('Power Spectral Density') 
for k = 1: length(SNR)
    
    for j = 1 : 20
        % generate white gaussian noise 
        noise = No/2 * randn(1,numberOfBits * Tb);
        time = 1: length(P_pulse);
        %  evaluate the Amplitude of the signal
        sigAmp = sqrt(power(10, (SNR(k)/10))* (2*No)/Tb);
%BPSK----------------------------------------------------------------------    
        % evaluate the BPSK signal before noise
        P_signal = P_pulse* sigAmp.* cos(Wc*time);

        % add noise to signal
        P_sig_Noise= P_signal + noise;

        % calculate H(t) of matched filter for bit 1
        P_Ht=0;
        for i = 1: Tb
            P_Ht(i) = sigAmp* cos(Wc*(Tb-i));
        end

        % evaluate the output of the matched filter
        P_matchOutput = conv(P_sig_Noise(1:40), P_Ht);
        P_decisionTaker = 0;
        P_count = 1;
        if (P_matchOutput(Tb) > 0)
              P_decisionTaker (P_count) = 1;
              P_count = P_count +1;
        else
              P_decisionTaker (P_count) = 0;
              P_count = P_count +1;
        end
        
        for i = 1:numberOfBits-1
            P_matchOutput = conv(P_sig_Noise(i*Tb + 1:(i+1)*Tb), P_Ht);
            if (P_matchOutput(Tb) > 0)
                P_decisionTaker (P_count) = 1;
                P_count = P_count +1;
            else
                P_decisionTaker (P_count) = 0;
                P_count = P_count +1;
            end
        end

        % get the total number of erroneous bits
        P_errorBits = 0;
        P_decisionTaker(P_decisionTaker == 0)= -1;
        for i = 1:numberOfBits
            if(P_decisionTaker(i) ~= BPSK_bits(i))
                P_errorBits = P_errorBits + 1;
            end
        end
        
        ErrorBPSK(j) = P_errorBits/numberOfBits;
%--------------------------------------------------------------------------
%BFSK----------------------------------------------------------------------
        % evaluate the BFSK signal
        OnePeriod = sigAmp*F_pulse.*cos(W1*time);
        newFpulse = F_pulse;
        %change the value of 1 to any number to change zero to 1 and return
        %1 to zero to evaluate the zero period
        newFpulse(newFpulse == 1)= 2;
        newFpulse(newFpulse == 0)= 1;
        newFpulse(newFpulse == 2)= 0;
        ZeroPeriod = sigAmp*newFpulse.*cos(W2*time);
        F_signal = OnePeriod + ZeroPeriod;
        % add noise to the signal 
        F_sig_Noise = F_signal + noise;
        % evaluate the h(t) for zero & one
        F_Ht_1=0;
        F_Ht_0=0;
        for i = 1: Tb
            F_Ht_1(i) = sigAmp* cos(W1*(Tb-i));
            F_Ht_0(i) = sigAmp* cos(W2*(Tb-i));
        end
        % get the output of the matched filter
        F_matchOutput_1 = conv(F_sig_Noise(1:40), F_Ht_1);
        F_matchOutput_0 = conv(F_sig_Noise(1:40), F_Ht_0);
        F_decisionTaker = 0;
        F_count = 1;
        if ((F_matchOutput_1(Tb) - F_matchOutput_0(Tb)) > 0)
              F_decisionTaker (F_count) = 1;
              F_count = F_count +1;
        else
              F_decisionTaker (F_count) = 0;
              F_count = F_count +1;
        end
        
        for i = 1:numberOfBits-1
            F_matchOutput_1 = conv(F_sig_Noise(i*Tb + 1:(i+1)*Tb), F_Ht_1);
            F_matchOutput_0 = conv(F_sig_Noise(i*Tb + 1:(i+1)*Tb), F_Ht_0);
            if ((F_matchOutput_1(Tb) - F_matchOutput_0(Tb)) > 0)
                F_decisionTaker (F_count) = 1;
                F_count = F_count +1;
            else
                F_decisionTaker (F_count) = 0;
                F_count = F_count +1;
            end
        end
        
        % get the total number of erroneous bits
        F_errorBits = 0;
       
        for i = 1:numberOfBits
            if(F_decisionTaker(i) ~= BFSK_bits(i))
                F_errorBits = F_errorBits + 1;
            end
        end
        
        ErrorBFSK(j) = F_errorBits/numberOfBits;
%--------------------------------------------------------------------------
    end
    AverageErrorBPSK(k) = sum(ErrorBPSK)/20;  
    AverageErrorBFSK(k) = sum(ErrorBFSK)/20;
    BPSK_TH_BE(k) = (1/2)*erfc(sigAmp*sqrt(Tb/(2*No)));
    BFSK_TH_BE(k) = (1/2)*erfc((sigAmp/2)*sqrt(Tb/No));
    
    
    if (SNR(k) == -2) 
        % draw generated bits, noise, modulated signal, modulatedsignal + noise
        % output of matched filter, the generated bits after decision taker 
        % for SNR = value in if and last generated noise for BPSK
        % draw pulses
        figure;
        plot (P_pulse);
        title('Input Pulses')
        xlabel('Time(sec)') 
        ylabel('Amplitude(v)')
        ylim([-1.5 1.5]);
        % draw noise
        figure;
        plot (noise);
        title('Noise')
        xlabel('Time(sec)') 
        ylabel('Amplitude(v)')
        % draw the signal after modulation
        figure;
        plot(P_signal);
        title('BPSK Modulated signal')
        xlabel('Time(sec)') 
        ylabel('Amplitude(v)')
        % draw the signal after modulation and adding the noise
        figure;
        plot(P_sig_Noise);
        title('BPSK Modulated signal + Noise')
        xlabel('Time(sec)') 
        ylabel('Amplitude(v)')
        % draw the generated bits after decision taker
        figure;
        P_newPulse = rectpulse(P_decisionTaker,Tb);
        plot(P_newPulse);
        ylim([-1.5 1.5]);
        title('Output of the matched filter')
        xlabel('Time(sec)') 
        ylabel('Amplitude(v)')

        % draw generated bits, noise, modulated signal, modulatedsignal + noise
        % output of matched filter, the generated bits after decision taker 
        % for SNR = value in if and last generated noise for BFSK
        % draw pulses
        figure;
        plot (F_pulse);
        ylim([-0.5 1.5]);
        title('Input Pulses')
        xlabel('Time(sec)') 
        ylabel('Amplitude(v)')
        % draw noise
        figure;
        plot (noise);
        title('Noise')
        xlabel('Time(sec)') 
        ylabel('Amplitude(v)')
        % draw signal after modulation
        figure;
        plot (F_signal);
        title('BFSK Modulated signal')
        xlabel('Time(sec)') 
        ylabel('Amplitude(v)')
        % draw the signal after modulation and adding the noise
        figure;
        plot(F_sig_Noise);
        title('BFSK Modulated signal + Noise')
        xlabel('Time(sec)') 
        ylabel('Amplitude(v)')
        % draw the generated bits after decision taker
        figure;
        F_newPulse = rectpulse(F_decisionTaker,Tb);
        plot(F_newPulse);
        ylim([-0.5 1.5]);
        title('Output of the matched filter')
        xlabel('Time(sec)') 
        ylabel('Amplitude(v)')  
    end
end


% draw the practical Average Bit Error Rate against Signal to Noise ratio 
figure 
plot (SNR,log(AverageErrorBPSK))
hold on;
plot (SNR,log(AverageErrorBFSK))
hold on; 
plot (SNR,log(BPSK_TH_BE)); 
hold on; 
plot (SNR,log(BFSK_TH_BE)); 
title('Average Bit Error Rate against Signal to Noise ratio')
xlabel('Signal to Noise ratio') 
ylabel('Bit Error Rate')
legend('BPSK', 'BFSK','BPSK Theoretical','BFSK Theoretical' )



