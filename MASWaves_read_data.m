
function [u,t,Tmax,L,x] = MASWaves_read_data(Filename,HeaderLines,fs,N,dx,x1,Direction)

% Read recorded data into matrix u
if strcmp(Direction,'forward')
    u = dlmread(Filename,'',HeaderLines,0);
elseif strcmp(Direction,'backward')
    u_temp = dlmread(Filename,'',HeaderLines,0); 
    for j = 1:N
        u(:,j) = u_temp(:,N-j+1);
    end
else
    disp('Error')
    u = [];
end

% Total recording time [s]
Tmax = length(u(:,1))/fs - 1/fs; 

% Time of individual recordings [s]
%T = linspace(0,Tmax,length(u(:,1))); 
t=1/fs:1/fs:length(u(:,1))/fs;
% Length of receiver spread [m]
L = (N-1)*dx; 

% Location of receivers, distance from seismic source [m]
x = (x1):dx:(L+x1);

end
