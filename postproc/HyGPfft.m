%%% script to compare the frequency spectrum (FFT) of original signal and
%%% HyGP model

function HyGPfft(R, tree_output, experiment, run, cur_gen)

% settings ---------------------------------------
size_title = 11;
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
T=0.0022
Fs=1/T
L=3000

% Original signal (input data set)
S = R(:,2);

% HyGP model
X=tree_output;

%Plot the HyGP model in the time domain
figure
plot(Z1,X)
% captions
line1_title = experiment;
line2_title = ['Sampled window of HyGP model'];
title({line1_title; line2_title}, 'FontSize',size_title,'Interpreter','none');
xlabel('Z1')
ylabel('X(Z1)')
xlim([1 11])

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
f = Fs*(0:(L/2))/L;  %if sampling frequency is Fs, fft can get up to Fs/2 
figure
plot(f,P1,'-k',f,P1_HyGP, '-r') 
line1_title = experiment;
line2_title =  ['Run ' num2str(run) ', Gen. ' num2str(cur_gen) ' - Building data set'];
line3_title = ['Single-Sided Amplitude Spectrum of HyGP model and original signal'];
title({line1_title; line2_title; line3_title}, 'FontSize',size_title,'Interpreter','none');
xlabel('f (Hz)')
ylabel('|P1(f)|')
legend('Original signal','HyGP model','Location','best')
