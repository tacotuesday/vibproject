%% Housekeeping
clear all
close all;

%% Creating H matrix
%load no_weight_data.mat;
load weight_data_3.mat;   % Comment/uncomment whichever data file is being analyzed.

[fft,freqH] = fft_func(data.x_sample,data.Nf);

F(1,:,:) = fft(1,:,:);
ZF = conj(F);

for p = 1:10
    for i = 1:7
        X(i,:,p) = fft(i+1,:,p);
        ZX = conj(X);
    end
end

xPowFX = 0.1.*sum(ZF.*X,3);
aPowFF = 0.1.*sum(ZF.*F,3);
xPowXF = 0.1.*sum(ZX.*F,3);
aPowXX = 0.1.*sum(ZX.*X,3);

H1 = xPowFX./aPowFF;
H2 = aPowXX./xPowXF;
H_est = (H2 + H1)./2;

coh = abs(H1./H2);
H_est(:,ceil(data.Nf+1):end,:) = []; % Truncates H matrix
H = mean(H_est, 3); % Average all trials to create H

%% Quadrature Code

n = size(H,2); % number of data points taken in H vector
freq_range = data.Nf; % Use the Nyquist frequency to determine the sampling frequency range (Hz)
freq_spacing = 1; % Frequency spacing (Hz) from the FFT transformation
freq = [0:(freq_range)/(n-1):freq_range]; % Create a frequency vector (Hz)

%% 1. Plotting Re[H] vs. omega to find undamped natural frequencies.
Hreal=real(H); % Real components of H

% plot routine
for k = 1:7
    figure(8)
    subplot(3,3,k);
    plot(freq,Hreal(k,:))
    title(['Real Component of Measured FRF H_' num2str(k) '_2'])
    grid on
    xlabel('frequency (Hz)'); 
    ylabel('Re[H(\omega)]'); % no units given for Re[H(\omega)]
end

figure(9)
plot(freq,Hreal)
title('Real Component of Measured FRF')
xlabel('frequency (Hz)'); 
ylabel('Re[H(\omega)]'); % no units given for Re[H(\omega)]
legend('H_1_2','H_2_2','H_3_2','H_4_2','H_5_2','H_6_2','H_7_2');
grid on

% Determine natural frequencies from FRF of each accelerometer
omega = [40.6, 79.6, 100, 150.1, 298, 451, 612];
% omega_1 = [147.8, 150.1, 152, 150.5, 152.3, 150.6, 150.5];
% omega_2 = [291.5, 302.7, 305.3, 295, 291.8, 291.8, 292.7];
% omega_3 = [450.8, 450.8, 451.0, 450.6, 456, 452.2, 451.7];
% omega(1) = mean(omega_1, 2); % Units: Hz
% omega(2) = mean(omega_2, 2);
% omega(3) = mean(omega_3, 2);
% fprintf('Measured natural frequencies (Hz):')
% disp(omega)



%% 2. finding damping ratios by the half-power point method
Himag = imag(H); % imaginary component of H
Hmag = sqrt(Hreal.^2.+Himag.^2); % magnitude of H

% the frequencies must be partitioned out into 3 sections to plot lines
spacing_bucket = [35, 60, 65, 85, 90, 110, 140, 160, 250, 350, 400, 500, 550, 675]; % Boundaries for each spacing bucket (Hz), determine this from natural frequencies above.
space1 = [35:60];
space2 = [65:85];
space3 = [90:110];
space4 = [140:160];
space5 = [250:350];
space6 = [400:500];
space7 = [550:675];
freq1 = [spacing_bucket(1):1/spacing_bucket(2):spacing_bucket(2)]; 
freq2 = [spacing_bucket(3):1/(spacing_bucket(4)-spacing_bucket(3)):spacing_bucket(4)];
freq3 = [spacing_bucket(5):1/(spacing_bucket(6)-spacing_bucket(5)):spacing_bucket(6)];
freq4 = [spacing_bucket(7):1/(spacing_bucket(8)-spacing_bucket(7)):spacing_bucket(8)];
freq5 = [spacing_bucket(9):1/(spacing_bucket(10)-spacing_bucket(9)):spacing_bucket(10)];
freq6 = [spacing_bucket(11):1/(spacing_bucket(12)-spacing_bucket(11)):spacing_bucket(12)];
freq7 = [spacing_bucket(13):1/(spacing_bucket(14)-spacing_bucket(13)):spacing_bucket(14)];

% for each section, take the half-power point amplitude
% for i = 1:7
%     for j = 1:3
%         hpp(i,j) = max(Hmag(i,space(j,1)))/sqrt(2);
%     end
% end

space1 = [35:60];
space2 = [65:85];
space3 = [90:110];
space4 = [140:160];
space5 = [250:350];
space6 = [400:500];
space7 = [550:675];

freqLength=0:20;

hpp(1,1) = max(Hmag(1,(space1)))/sqrt(2);
hpp(1,2) = max(Hmag(1,(space2)))/sqrt(2);
hpp(1,3) = max(Hmag(1,(space3)))/sqrt(2);
hpp(1,4) = max(Hmag(1,(space4)))/sqrt(2);
hpp(1,5) = max(Hmag(1,(space5)))/sqrt(2);
hpp(1,6) = max(Hmag(1,(space6)))/sqrt(2);
hpp(1,7) = max(Hmag(1,(space7)))/sqrt(2);

hpp(2,1) = max(Hmag(2,(space1)))/sqrt(2);
hpp(2,2) = max(Hmag(2,(space2)))/sqrt(2);
hpp(2,3) = max(Hmag(2,(space3)))/sqrt(2);
hpp(2,4) = max(Hmag(2,(space4)))/sqrt(2);
hpp(2,5) = max(Hmag(2,(space5)))/sqrt(2);
hpp(2,6) = max(Hmag(2,(space6)))/sqrt(2);
hpp(2,7) = max(Hmag(2,(space7)))/sqrt(2);

hpp(3,1) = max(Hmag(3,(space1)))/sqrt(2);
hpp(3,2) = max(Hmag(3,(space2)))/sqrt(2);
hpp(3,3) = max(Hmag(3,(space3)))/sqrt(2);
hpp(3,4) = max(Hmag(3,(space4)))/sqrt(2);
hpp(3,5) = max(Hmag(3,(space5)))/sqrt(2);
hpp(3,6) = max(Hmag(3,(space6)))/sqrt(2);
hpp(3,7) = max(Hmag(3,(space7)))/sqrt(2);

hpp(4,1) = max(Hmag(4,(space1)))/sqrt(2);
hpp(4,2) = max(Hmag(4,(space2)))/sqrt(2);
hpp(4,3) = max(Hmag(4,(space3)))/sqrt(2);
hpp(4,4) = max(Hmag(4,(space4)))/sqrt(2);
hpp(4,5) = max(Hmag(4,(space5)))/sqrt(2);
hpp(4,6) = max(Hmag(4,(space6)))/sqrt(2);
hpp(4,7) = max(Hmag(4,(space7)))/sqrt(2);

hpp(5,1) = max(Hmag(5,(space1)))/sqrt(2);
hpp(5,2) = max(Hmag(5,(space2)))/sqrt(2);
hpp(5,3) = max(Hmag(5,(space3)))/sqrt(2);
hpp(5,4) = max(Hmag(5,(space4)))/sqrt(2);
hpp(5,5) = max(Hmag(5,(space5)))/sqrt(2);
hpp(5,6) = max(Hmag(5,(space6)))/sqrt(2);
hpp(5,7) = max(Hmag(5,(space7)))/sqrt(2);

hpp(6,1) = max(Hmag(6,(space1)))/sqrt(2);
hpp(6,2) = max(Hmag(6,(space2)))/sqrt(2);
hpp(6,3) = max(Hmag(6,(space3)))/sqrt(2);
hpp(6,4) = max(Hmag(6,(space4)))/sqrt(2);
hpp(6,5) = max(Hmag(6,(space5)))/sqrt(2);
hpp(6,6) = max(Hmag(6,(space6)))/sqrt(2);
hpp(6,7) = max(Hmag(6,(space7)))/sqrt(2);

hpp(7,1) = max(Hmag(7,(space1)))/sqrt(2);
hpp(7,2) = max(Hmag(7,(space2)))/sqrt(2);
hpp(7,3) = max(Hmag(7,(space3)))/sqrt(2);
hpp(7,4) = max(Hmag(7,(space4)))/sqrt(2);
hpp(7,5) = max(Hmag(7,(space5)))/sqrt(2);
hpp(7,6) = max(Hmag(7,(space6)))/sqrt(2);
hpp(7,7) = max(Hmag(7,(space7)))/sqrt(2);


for p = 1:7
    figure(p) % |H| vs frequency
    semilogy(freq,Hmag(p,:)); 
    hold on;
    plot(freq1,hpp(p,1)*ones(1,length(freq1)));
    plot(freq2,hpp(p,2)*ones(1,length(freq2)));
    plot(freq3,hpp(p,3)*ones(1,length(freq3)));
    plot(freq4,hpp(p,5)*ones(1,length(freq4)));
    plot(freq5,hpp(p,5)*ones(1,length(freq5)));
    plot(freq6,hpp(p,6)*ones(1,length(freq6)));
    plot(freq7,hpp(p,7)*ones(1,length(freq7)));
    % plot(freq7,hpp(p,3)*ones(1,length(freq7)));
    title('Magnutude of measured FRFs of a force applied to 2nd mass')
    xlabel('frequency (Hz)'); ylabel('H(\omega)'); % no units given for Re[H(\omega)]
    legend(['H_' num2str(p) '_2'],['hpp' num2str(p) '_1'],['hpp' num2str(p) '_2'],...
        ['hpp' num2str(p) '_3']);
    grid on; 
end

omegazeta(1,1,:) = [36.52, 46.5];
omegazeta(1,2,:) = [72.25, 80.6];
omegazeta(1,3,:) = [96.8, 107.3];
omegazeta(1,4,:) = [NaN, NaN];
omegazeta(1,5,:) = [288.7, 307];
omegazeta(1,6,:) = [436, 457.6];
omegazeta(1,7,:) = [596.9, 624.3];
omegazeta(2,1,:) = [38.7, 45.12];
omegazeta(2,2,:) = [77.5, 81.75];
omegazeta(2,3,:) = [96.1, 105.8];
omegazeta(2,4,:) = [144.3, 161.1];
omegazeta(2,5,:) = [295.8, 319.8];
omegazeta(2,6,:) = [442.9, 463.6];
omegazeta(2,7,:) = [604.1, 627];
omegazeta(3,1,:) = [38.02, 42.93];
omegazeta(3,2,:) = [77.35, 82.25];
omegazeta(3,3,:) = [92.8, 103.8];
omegazeta(3,4,:) = [NaN, NaN];
omegazeta(3,5,:) = [295.9, 318.1];
omegazeta(3,6,:) = [441.3, 460.3];
omegazeta(3,7,:) = [653.9, 685.3];
omegazeta(4,1,:) = [36.38, 43.95];
omegazeta(4,2,:) = [78.05, 83.2];
omegazeta(4,3,:) = [NaN, 102.2];
omegazeta(4,4,:) = [143.5, 156.8];
omegazeta(4,5,:) = [285.5, 299.9];
omegazeta(4,6,:) = [443.5, 464.4];
omegazeta(4,7,:) = [601.8, 624.7];
omegazeta(5,1,:) = [39.65, 44.68];
omegazeta(5,2,:) = [78.35, NaN];
omegazeta(5,3,:) = [NaN, 102.3];
omegazeta(5,4,:) = [NaN, NaN];
omegazeta(5,5,:) = [278.6, 328];
omegazeta(5,6,:) = [438.9, 463.1];
omegazeta(5,7,:) = [608.3, 666.8];
omegazeta(6,1,:) = [38.32, 43.4];
omegazeta(6,2,:) = [75.55, 81];
omegazeta(6,3,:) = [95.9, 105.7];
omegazeta(6,4,:) = [NaN, NaN];
omegazeta(6,5,:) = [287.5, 306.1];
omegazeta(6,6,:) = [443.1, 483.8];
omegazeta(6,7,:) = [649.8, NaN];
omegazeta(7,1,:) = [39.63, 44.4];
omegazeta(7,2,:) = [78.1, 83.25];
omegazeta(7,3,:) = [NaN, NaN];
omegazeta(7,4,:) = [NaN, NaN];
omegazeta(7,5,:) = [279.6, 304];
omegazeta(7,6,:) = [440.4, 461.7];
omegazeta(7,7,:) = [605.2, 635.8];

for i = 1:7 		% Calculating possible zeta values
    for j = 1:7
        zetar(i,j,1) = sqrt(1/2 + sqrt(1 - (((omegazeta(i,j,1)^2 - omegazeta(i,j,2)^2)/omega(j)^2)^2)/4)/2);
        zetar(i,j,2) = sqrt(1/2 - sqrt(1 - (((omegazeta(i,j,1)^2 - omegazeta(i,j,2)^2)/omega(j)^2)^2)/4)/2);
    end
end

for i = 1:7 	% rewriting averaging function into a general case to eliminate bad values
    for j = 1:7
        for k = 1:2 	% positive and negative value only
            if zetar(i,j,k)<0.7 % assuming smaller value is good
                zetaverage(i,j) = zetar(i,j,k);
            else
                zetaverage(i,j) = 0;
            end
        end
    end
end

for i = 1:7 % averaging zeta values for final result
    zeta(i) = mean(nonzeros(zetaverage(:,i)));
end

% fprintf('Measured damping ratios:')
% disp(zeta)

%% 3. Finding the mode shapes by plotting Im[H] vs. omega, and modal masses

% plot routine
figure(11)   % Im[H] vs frequency
plot(freq, Himag);
for i = 1:3
    line([omega(i) omega(i)], [-3 3]); % vertical lines at each natural frequency
end
title('Imag component measured FRFs of a force applied to 2nd mass')
xlabel('frequency (rad/s)'); ylabel('Im[H(\omega)]'); % no units given for Im[H(\omega)]
legend('H_1_2', 'H_2_2', 'H_3_2', 'H_4_2', 'H_5_2','H_6_2','H_7_2'); grid;


Himag1 = [.09309 .8998 -.1101 -.804 -2.258 .2747 2.667];
Himag2 = [.325 1.213 -1.032 -.7456 .2156 .6487 -.9926];
Himag3 = [-.7418 -1.509 -1.794 -.7735 -2.442 -1.799 -2.917];
Himag4 = [.7001 .2921 -.2832 .6658 -.1665 2.192 -.4606];
Himag5 = [-1.292 -.03444 .522 -.5505 2.043 1.855 3.191];
Himag6 = [.9439 -1.049 1.161 -.2492 -.1494 .1181 1.21];
Himag7 = [-.3993 .7091 .06277 .4313 -2.95 .1746 2.341];
Hmode = [Himag1; Himag2; Himag3; Himag4; Himag5; Himag6; Himag7];

% Ar, the residue matrix can be calculated since zeta, omega, and FRF magnitudes are known
% Ar is only calculated for measured modes, and is formatted accordingly
for i = 1:7
    for r = 1:7
        Ar(i,r) = -2*zeta(r)*(omega(r)*2*pi)^2.*Hmode(i,r);
    end
end

u_measured(1,1:7) = 1; % scaling assumption

for r = 1:7
    m(r) = 1/Ar(1,r); % modal masses can be calculated here assuming response x1=1
    for i = 1:7
        u_measured(i,r) = Ar(i,r).*m(r)./u_measured(1,r); % The rest of the u matrix can also be calculated
    end
end

% fprintf('Measured modal masses (assuming x1=1):')
disp(m)


%% PLOT MODE SHAPES: Original loop would have been great, but Richards wants
% one figure per mode shape.
% Plot Mode 1
figure(22)
plot(1:7,[0 0 0 0 0 0 0]);
hold on
plot(1:7, u_measured(:,1))
title('Mode 1 Shape (Modified)')
xlabel('Location'); ylabel('Relative Magnitude'); 
axis([1 7 -17 17]);
legend('Mode 1','Original Shape');

% Plot Mode 2
figure(23)
plot(1:7,[0 0 0 0 0 0 0]);
hold on
plot(1:7, u_measured(:,2))
title('Mode 2 Shape (Modified)')
xlabel('Location'); ylabel('Relative Magnitude'); 
axis([1 7 -17 17]);
legend('Mode 2','Original Shape');

% Plot Mode 3
figure(24)
plot(1:7,[0 0 0 0 0 0 0]);
hold on
plot(1:7, u_measured(:,3))
title('Mode 3 Shape (Modified)')
xlabel('Location'); ylabel('Relative Magnitude'); 
axis([1 7 -17 17]);
legend('Mode 3','Original Shape');
%% Synthesized Frequency Response Functions
% Much of this code is adapated from the
% superposition_freq_domain_altered.m class example

Ar(:,:,2) = Ar; % properly setting the residue matrix
Ar(:,:,1) = 0;
% FRF matrix synthesis
for cc1 = 1:length(freqH)
    for i = 1:7
        for r = 1:7
            for j = 1:7
               A(i,r,j) = u_measured(i,r)*u_measured(j,r).'/m(r);
               H_h(i,r,j,cc1) = A(i,r,j)./(omega(r)^2 - (2*pi*freqH(cc1))^2 + 2*1i*zetar(r)*omega(r)*2*pi*freqH(cc1));
            end
        end
    end
end

HT = squeeze(sum(H_h,3)); % superimpose modes => total "measured" FRF matrix

for i = 1:7
    for j = 1:7
        for r = 1:7
            H_plot  = squeeze(H_h(i,r,:,:));
            HT_synth_helper1 = squeeze(HT(i,j,:));
            HT_synth_helper2 = transpose(HT_synth_helper1);
            HT_synth (i,:,r) = HT_synth_helper2;
        end
    end
end

% delete this plot routine later
for pm = 1:7 % plotting synthesized FRFs with measured data - magnitude and phase
    figure()
    subplot(2,1,1);
    plot(abs(HT_synth(pm,:,2))); hold on
    title(['Sensor Location ' num2str(pm) ' Responses'])
    xlabel('Frequency (Hz)'); ylabel('|H|');
    subplot(2,1,2);
    plot(freqH, atan2(-imag(HT_synth(pm,:,2)), real(HT_synth(pm,:,2)))*180/pi);
    xlabel('Frequency (Hz)'); ylabel('Phase (deg)');
    hold off
end

%% Comparing Analytical and Measured methods - Creating the output plots

% Nathan
% for i = 1:3 % mode shape plots
%     figure();
%     plot(1:5, u_measured(:,i)); hold on % measured mode shape
%     plot(1:4/(length(u_analytical)-1):5, u_analytical(i,:)); % analytical mode shape
%     % title(['Sensor Location ' num2str(i) ' Mode Shape'])
%     xlabel('Location'); ylabel('Relative Intensity'); set(gca,'XTick',0:5);
%     ylim([-1 1]);
%     legend('Measured', 'Analytical'); grid; hold off;
% end
% disp('mode shape plots have been created');

%% Creating Measured FRF and Coherence plots
for pl = 1:7 % plotting coherence and FRFs
    figure()
    subplot(2,1,1);
    plot(freqH,abs(X(pl,:,1))); hold on
    title(['Sensor Location ' num2str(pl) ' Response (Unmodified)'])
    xlabel('Frequency (Hz)'); ylabel('|H|');
    subplot(2,1,2);
    plot(freqH,coh(pl,:));
    xlabel('Frequency (Hz)'); ylabel('Coherence');
    hold off
end


%THIS COULD BE USEFUL FOR REPORT BUT WAY TOO MANY PLOTS SO I COMMENTED OUT
% for pm = 1:5 % plotting synthesized FRFs with measured data - magnitude and phase
%     figure()
%     subplot(2,1,1);
%     plot(freqH,abs(X(pm,:,1))); hold on
%     plot(freqH,abs(HT_synth(pm,:,2)));
%     title(['Sensor Location ' num2str(pm) ' Responses'])
%     xlabel('Frequency (Hz)'); ylabel('|H|');
%     legend('Measured', ' Synthesized');
%     subplot(2,1,2);
%     plot(freqH, atan2(-imag(HT_synth(pm,:,2)), real(HT_synth(pm,:,2)))*180/pi);
%     xlabel('Frequency (Hz)'); ylabel('Phase (deg)');
%     hold off
% end
% 
% for i = 1:7 % plotting synthesized FRFs with measured data - magnitude and phase
%     for k = 1:3
%         figure()
%         subplot(2,1,1);
%         plot(freqH,abs(HT_synth(i,:,2)));
%         title(['Force Location ' num2str(k) ' Sensor Location ' num2str(i) ' Responses'])
%         xlabel('Frequency (Hz)'); ylabel('|H|');
%         subplot(2,1,2);
%         plot(freqH, atan2(-imag(HT_synth(i,:,k)), real(HT_synth(i,:,k)))*180/pi);
%         xlabel('Frequency (Hz)'); ylabel('Phase (deg)');
%         hold off
%     end
% end
