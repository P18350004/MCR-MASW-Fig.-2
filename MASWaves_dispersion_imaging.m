
function [f,c,A] = MASWaves_dispersion_imaging(u,N,x,fs,cT_min,cT_max,delta_cT)

% Converting measuring frequency from Hz to rad/sec
omega_fs = 2*pi*fs; 

% Number of samples in each trace
Lu = length(u(:,1));

% Empty matrices with Lu lines and n columns
U = zeros(Lu,N);
P = zeros(Lu,N);
Unorm = zeros(Lu,N);

% Apply discrete Fourier transform to the time axis of u
% (The MATLAB function fft( ) is used)
for j = 1:N
    U(:,j) = fft(u(:,j));
end

% Number of samples in each transformed trace
LU = length(U(:,1));

% Normalize U in offset and frequency domains
% Compute the phase spectrum of U
i = sqrt(-1);
for j = 1:N
    for k = 1:LU
        Unorm(k,j) = U(k,j)/abs(U(k,j));
    end
    P(:,j) = exp(i.*-angle(U(:,j)));
end

% Frequency range for U
% Lomega = LU
omega = (0:LU-1)*(omega_fs/LU); 

% Compute the slant-stack (summed) amplitude corresponding to each set of
% omega and cT, A(omega,cT).
cT = cT_min:delta_cT:cT_max;
LcT = length(cT);

% Empty matrices with LU lines and LcT columns
c = zeros(LU,LcT);
f = zeros(LU,LcT);
A = zeros(LU,LcT); 

for j = 1:LU % Frequency component j
    for k = 1:LcT % Testing Rayleigh wave phase velocity k
        % Frequency (in [Hz]) corresponding to angular frequency omega
        f(j,k) = omega(j)/(2*pi);
        % Testing phase velocity [m/s]
        c(j,k) = cT(k);
        
        delta = omega(j)/cT(k);
       
        temp = 0;
        for l = 1:N
            temp = temp + exp(-i*delta*x(l))*P(j,l);     
        end
        % Compute absolute value and normalize with respect to number of 
        % receivers
        A(j,k) = abs(temp)/N; 
    end
end
end