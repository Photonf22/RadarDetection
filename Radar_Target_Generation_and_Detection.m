clear all
clc;
%dbstop in RadarDetectionFinalProject.m at 156
%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
% velocity can be any value in the range of -70 to + 70 m/s.
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%speed of light = 3e8
%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains contant
v_init= -20; %m/s initial velocity of the target
R_init= 110; % initial range of the target
Vmax= 100;
Rmax= 200;
%% FMCW Waveform Generation

% *%TODO* :
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% Bandwidth(Bsweep)=speedoflight/(2∗rangeResolution)
% chirp using the requirements above.
c=  3*10^8;         % speed of light
range_resolution = 1;
Tchirp=5.5*2*Rmax/c; % chirp time or TChirp
%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq
% TODO: Calculate the wavelength
wav= c/fc;            % wavelength of carrier signal
% B(sweep) of chirp or bandwidth
B_sweep= c /(2*range_resolution);

Slope =B_sweep/Tchirp ;   %Slope (slope) of the FMCW
disp(Slope)
disp(B_sweep)
disp(Tchirp)
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples


% Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));

%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 
for i=1:length(t)        
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity.
    % t(i) - t(i-1); trip time/ difference between current time and
    % previous time would give us the trip time thus so far
    
    % calculate current distance and add it to initial distance for every
    % loop to get the total distance thus so far per loop which then put
    % into r_t(i)   
    if i == 1                      % if first loop then we are @ init dist
        r_t(i) = R_init;
    else
        r_t(i) = r_t(i-1) + (v_init*(t(i) - t(i-1)));
    end
    % td is the time delay in which R is the distance thus traveled so far
    % which is 'r_t' which represents the distance at this current moment
    % overall from the target's initial distance
    % td= 2*R/c; where R is the distance traveled so far
    td(i) =2*r_t(i) / c;
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal.
    Tx(i) = cos(2*pi*((fc*t(i)) + (Slope*(t(i)^2))/2));
    Rx(i) = cos(2*pi*(fc*(t(i)- td(i)) + (Slope*((t(i) - td(i))^2)/2)));
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i).*Rx(i);
end

%% RANGE MEASUREMENT
 % *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.

sig  = reshape(Mix, [Nr, Nd]);
 % *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
 % *%TODO* :
% Take the absolute value of FFT output
sig_fft = fft(sig,Nr);
% normalizing the fft signal
sig_fft_norm = (sig_fft)/max(sig_fft);
% taking absolute value of the signal
sig_fft_sub = abs(sig_fft_norm);
 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
sig_fft_half = sig_fft_sub(1:Nr/2 - 1);
%plotting the range
figure ('Name','Range from First FFT')
subplot(2,1,1)

 % *%TODO* :
 % plot FFT output
 % TODO : Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
% Plotting
plot(sig_fft_half);
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
axis ([0 200 0 1]);

%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

CFAR = zeros(size(RDM));
%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);
%% CFAR implementation

%Slide Window through the complete Range Doppler Map

% *%TODO* :
%Select the number of Training Cells in both the dimensions.
Tr=11;
Td=9;
Gr=5;
Gd=3;
Grid_size= (2*Tr+2*Gr+1)*(2*Td+2*Gd+1);
T =  (2*Tr+2*Gr+1)*(2*Td+2*Gd+1) - (2*Gr+1)*(2*Gr+1);
% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
G = (2*Gr+1)*(2*Gd+1) - 1;
% *%TODO* :
% offset the threshold by SNR value in dB
offset=12;
% *%TODO* :
%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(1,1);

% *%TODO* :
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.
% Slide window across the signal length
max_T =1;
[m,n] = size(RDM) 
for i = Td+Gd+1:m-(Gd+Td)
    for j = Tr+Gr+1:n-(Gr+Tr)
        noise_level= 0;
            for p = i - (Td+Gd):i+Td+Gd
                for q = j-(Tr+Gr):j+Tr+Gr
                    if(abs(i-p)>Gd && abs(j-q)>Gr)
                        noise_level= noise_level + db2pow(RDM(p,q)); 
                    end
                end   
            end
        threshold= pow2db(noise_level/(2*(Td+Gd+1)*2*(Tr+Gr+1)-(Gr*Gd)-1));
        threshold= threshold + offset;
% *%TODO* :
% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 
        CUT=RDM(i,j);        
        if(CUT < threshold)
            CFAR(i,j) = 0;
        else 
            CFAR(i,j) = max_T;
        end
    end
end
   % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
   % CFAR

% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure,surf(doppler_axis,range_axis,CFAR);
colorbar;
