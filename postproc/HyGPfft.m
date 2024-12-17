%%% script to compare the frequency spectrum (FFT) of original signal and
%%% HyGP model

function HyGPfft(R, tree_output, experiment, run, cur_gen, dataset_type)

% if filtering is used, update passband central frequency manually (line 94)! 

% settings ---------------------------------------
size_title = 14;
size_axes = 12;
pi=3.141592653589793e+00;
sep = filesep;  % set separator 
%-------------------------------------------------------

% R is the input data set, so:
% first column= independent variable
% second column= value of the original signal 
R(1:10,:)  % just to check that the right data are used

% Specify the parameters of a signal with a sampling frequency of 1 kHz and a signal duration of 1.5 seconds.
%Fs = 454.5454; %to have sampling window equal to building data set (1-11)      % Sampling frequency                    
%T = 1/Fs;             % Sampling period       
%L = 5000;             % Length of signal - size of building data set?
%Z1 = ones(1,L)+(0:L-1)*T;        % Time vector on which sampling is made (to be replaced by input data set)
Z1=R(:,1);
disp('R(2,1)')
R(2,1)
disp('R(1,1)')
R(1,1)
T=R(2,1)-R(1,1) % trying to capture the delta between time steps - assuming that sampling is UNIFORM    %T=0.0022
Fs=1/T
L=size(R,1)                 %3000

% Original signal (input data set)
S = R(:,2);

% HyGP model
X=tree_output;

%Plot the HyGP model in the time domain
figure
plot(Z1,X)
% captions
line1_title = experiment;
line2_title = ['Sampled window of HyGP model - ' dataset_type ' data set'];
title({line1_title; line2_title}, 'FontSize',size_title,'Interpreter','none');
xlabel('Z1','FontSize',size_axes)
ylabel('X(Z1)','FontSize',size_axes)
%xlim([1 11])

% Compute the Fourier transform of the HyGP model.
Y_HyGP = fft(X);
% Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
P2_HyGP = abs(Y_HyGP/L);   % two-sided spectrum
P1_HyGP = P2_HyGP(1:L/2+1); % one-sided spectrum, half magnitude of frequency component (only the positive)
P1_HyGP(2:end-1) = 2*P1_HyGP(2:end-1); % % one-sided spectrum, correct magnitude 

% Fourier transform of the original signal
Y = fft(S);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

% Define the frequency domain f and plot the single-sided amplitude spectrum P1
f_P1 = Fs*(0:(L/2))/L;  %if sampling frequency is Fs, fft can get up to Fs/2 
figure
plot(f_P1,P1,'-k',f_P1,P1_HyGP, '-r') 
line1_title = experiment;
line2_title =  ['Run ' num2str(run) ', Gen. ' num2str(cur_gen) ' - ' dataset_type ' data set'];
line3_title = ['Single-Sided Amplitude Spectrum of HyGP model and original signal'];
title({line1_title; line2_title; line3_title}, 'FontSize',size_title,'Interpreter','none');
xlabel('f (Hz)')
ylabel('|P1(f)|')
legend('Original signal','HyGP model','Location','best')

% compute energy spectral density 
en_density = P2.^2;     % first square...
en_density_HyGP = P2_HyGP.^2;     % first square...
p_P1 = en_density(1:L/2+1);  
p_HyGP_P1 = en_density_HyGP(1:L/2+1);
p_P1(2:end-1) = 2*p_P1(2:end-1);  % then sum! For 15Hz component power = (0.5*0.5)*2 = 0.5  (https://en.wikipedia.org/wiki/Spectral_density)
p_HyGP_P1(2:end-1) = 2*p_HyGP_P1(2:end-1); 
figure()
plot(f_P1,p_P1,'-k',f_P1,p_HyGP_P1,'-r')
line1_title = experiment;
line2_title =  ['Run ' num2str(run) ', Gen. ' num2str(cur_gen) ' - ' dataset_type ' data set'];
line3_title = ['One-sided energy spectral density of HyGP model and original signal'];
title({line1_title; line2_title; line3_title}, 'FontSize',size_title,'Interpreter','none');
xlabel('f (Hz)')
ylabel('|P1(f)|^2')
legend('Original signal','HyGP model','Location','best')

% filter the signal power with a rectangular band passing filter
%f_central = 1787 %St=0.1 %3574 % St=0.2  
%f_amplitude= 0.05
%df = f_amplitude*f_central
%mask=(f_P1>=f_central-df) & (f_P1<=f_central+df);
mask=(f_P1>=180) & (f_P1<=17000);   % 17000 Hz corresponds to St=1.0
filtered_f_P1=f_P1(mask)
filtered_p_P1=p_P1(mask)
filtered_p_HyGP_P1=p_HyGP_P1(mask)
figure()
plot(filtered_f_P1,filtered_p_P1,'-k',filtered_f_P1,filtered_p_HyGP_P1,'-r')
line1_title = experiment;
line2_title =  ['Run ' num2str(run) ', Gen. ' num2str(cur_gen) ' - ' dataset_type ' data set'];
line3_title = ['Filtered one-sided signal energy of HyGP model and original signal (180-17000 Hz fr. range)'];
title({line1_title; line2_title; line3_title}, 'FontSize',size_title,'Interpreter','none');
xlabel('f (Hz)')
ylabel('|P1(f)|^2')
legend('Original signal','HyGP model','Location','best')

% integrate the filtered signal power to get energy
power= trapz(filtered_f_P1,filtered_p_P1)
power_HyGP= trapz(filtered_f_P1,filtered_p_HyGP_P1)
figure()
x=1;
b=bar(x,[power power_HyGP])
b(1).FaceColor = [0 0 0];
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')
b(2).FaceColor = [1 0 0];
xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom')
%ylim([0 max(power)*1.5])
line1_title = experiment;
line2_title =  ['Run ' num2str(run) ', Gen. ' num2str(cur_gen) ' - ' dataset_type ' data set'];
line3_title = ['Power of HyGP model and original signal after filtering'];
title({line1_title; line2_title; line3_title}, 'FontSize',size_title,'Interpreter','none');
legend('Original signal','HyGP model','Location','best')

% % % additional metrics based on relative error on maxima (Sergey 13/11/23)
% % %------------------------------------------------------------------------
% % % only for Towne's SPOD analysis for article stemming from ICAS paper (presented September 2022)
% % % set central frequency and frequency range
% % St_ref_frequency = input("Reference Strouhal frequency (Hz)? ")
% % freq_range_1=[0.9 1.1]*St_ref_frequency
% % freq_range_2=[0.8 1.2]*St_ref_frequency
% % 
% % % Find maxima in frequency range
% % % f_P1 entire frequency range
% % % P1  original signal FFT
% % % P1_HyGP   GP model FFT
% % % frequency range 1
% % mask1=(f_P1>=freq_range_1(1)) & (f_P1<=freq_range_1(2)); 
% % filtered_f_P1_1=f_P1(mask1)
% % filtered_P1_1=P1(mask1)
% % filtered_P1_HyGP_1=P1_HyGP(mask1)
% % max_P1_range1=max(filtered_P1_1)
% % max_P1_HyGP_range1=max(filtered_P1_HyGP_1)
% % % frequency range 2
% % mask2=(f_P1>=freq_range_2(1)) & (f_P1<=freq_range_2(2)); 
% % filtered_f_P1_2=f_P1(mask2)
% % filtered_P1_2=P1(mask2)
% % filtered_P1_HyGP_2=P1_HyGP(mask2)
% % max_P1_range2=max(filtered_P1_2)
% % max_P1_HyGP_range2=max(filtered_P1_HyGP_2)
% % 
% % % compute relative errors
% % err1=abs(max_P1_HyGP_range1-max_P1_range1)/abs(max_P1_range1)
% % err2=abs(max_P1_HyGP_range2-max_P1_range2)/abs(max_P1_range2)
% % 
% % %plots
% % figure()
% % plot(filtered_f_P1_1,filtered_P1_1,'-k',filtered_f_P1_1,filtered_P1_HyGP_1,'-r')
% % line1_title = experiment;
% % line2_title =  ['Run ' num2str(run) ', Gen. ' num2str(cur_gen) ' - ' dataset_type ' data set'];
% % line3_title = ['Relative error on [0.9 1.1] f_St = [' num2str(freq_range_1(1)) ' ' num2str(freq_range_1(2)) '] -> ' num2str(err1)];
% % title({line1_title; line2_title; line3_title}, 'FontSize',size_title,'Interpreter','none');
% % xlabel('f (Hz)')
% % ylabel('|P1(f)|')
% % legend('Original signal','HyGP model','Location','best')
% % 
% % figure()
% % plot(filtered_f_P1_2,filtered_P1_2,'-k',filtered_f_P1_2,filtered_P1_HyGP_2,'-r')
% % line1_title = experiment;
% % line2_title =  ['Run ' num2str(run) ', Gen. ' num2str(cur_gen) ' - ' dataset_type ' data set'];
% % line3_title = ['Relative error on [0.8 1.2] f_St = [' num2str(freq_range_2(1)) ' ' num2str(freq_range_2(2)) '] -> ' num2str(err2)];
% % title({line1_title; line2_title; line3_title}, 'FontSize',size_title,'Interpreter','none');
% % xlabel('f (Hz)')
% % ylabel('|P1(f)|')
% % legend('Original signal','HyGP model','Location','best')





