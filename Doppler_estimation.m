% Doppler Velocity Calculation
c = 3*10^8;         %speed of light
frequency = 77e9;   %frequency in Hz

% TODO: Calculate the wavelength
wav= c/frequency;

% TODO: Define the doppler shifts in Hz using the information from above 
fd_s1= 3e3;
fd_s2= -4.5e3; 
fd_s3= 11e3; 
fd_s4= -3e3;

% TODO: Calculate the velocity of the targets  fd = 2*vr/lambda
vr1= (fd_s1*wav)/2;
vr2= (fd_s2*wav)/2;
vr3= (fd_s3*wav)/2;
vr4= (fd_s4*wav)/2;
% TODO: Display results
disp(vr1);
disp(vr2);
disp(vr4);
disp(vr1);

% 200 m range ahead of the ego (radar)
% velocity of the ego vehicle is 5 m/s
% Then in next 5 seconds which of the targets would be closest to the ego 
% vehicle. Let%s give targets a tag [A B C D]
fd_radar= (2*5)/c; 
x_a= (200)+(5*(fd_s1 - fd_radar));
x_b= (200)+(5*(fd_s2 - fd_radar));
x_c= (200)+(5*(fd_s3 - fd_radar));
x_d= (200)+(5*(fd_s4 - fd_radar));

distance_total=[x_a,x_b,x_c,x_d];
disp(distance_total)