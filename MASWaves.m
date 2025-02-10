clear all
close all
clc
tic
%% DATA INPUT 
Filename = 'Shot4_0.3s.dat'; % wavefield file is given in respective folder
HeaderLines = 7;
fs = 4000; % Hz   Sampling Frequency 
N = 24;  % Number of Geophones
x1 = 5; % m  Offset Distance
dx = 1; % m   Geophone Spacing
Direction = 'forward';

[u,t,Tmax,L,x] = MASWaves_read_data(Filename,HeaderLines,fs,N,dx,x1,Direction);

%% PLOT DATA 
du = 1/25;
FigWidth = 20; % cm
FigHeight = 10; % cm
FigFontSize = 8; % pt

figure
MASWaves_plot_data(u,N,dx,x1,L,t,Tmax,du,FigWidth,FigHeight,FigFontSize)

%% DISPERSION   cT is the Rayliegh wave velocity
cT_min = 0; % m/s
cT_max = 3000; % m/s
delta_cT = 50; % m/s

[f,c,A] = MASWaves_dispersion_imaging(u,N,x,fs,cT_min,cT_max,delta_cT);

%% VELOCITY SPECTRA 2D
resolution = 100; 
fmin = 0; % Hz
fmax = 50; % Hz
FigWidth = 7; % cm
FigHeight = 7; % cm
FigFontSize = 8; % pt
figure
[fplot,cplot,Aplot] = MASWaves_plot_dispersion_image_2D(f,c,A,fmin,fmax,...
    resolution,FigWidth,FigHeight,FigFontSize);

%% VELOCITY SPECTRA 3D
fmin = 1; % Hz
fmax = 50; % Hz
FigWidth = 10; % cm
FigHeight = 10; % cm
FigFontSize = 8; % pt
figure
[fplot,cplot,Aplot] = MASWaves_plot_dispersion_image_3D(f,c,A,fmin,fmax,...
    FigWidth,FigHeight,FigFontSize);

%% SELECT POINTS 
f_receivers = 4.5; % Hz
 select = 'numbers'; % enter "1:14" after running code infront of "Fundamental mode dispersion curve:" in command window
% select = 'both';
up_low_boundary = 'no'; 
p = 95; % Percentage
[f_curve0,c_curve0,lambda_curve0,...
    f_curve0_up,c_curve0_up,lambda_curve0_up,...
    f_curve0_low,c_curve0_low,lambda_curve0_low] = ...
    MASWaves_extract_dispersion_curve(f,c,A,fmin,fmax,f_receivers,...
    select,up_low_boundary,p);

%% DISPERSION CURVE GENERATION
FigWidth = 9; % cm
FigHeight = 6; % cm
FigFontSize = 8; % pt
type = 'f_c';
up_low_boundary = 'yes';
figure
MASWaves_plot_dispersion_curve(f_curve0,c_curve0,lambda_curve0,...
     f_curve0_up,c_curve0_up,lambda_curve0_up,f_curve0_low,c_curve0_low,...
     lambda_curve0_low,type,up_low_boundary,FigWidth,FigHeight,FigFontSize)

FigWidth = 7; % cm
FigHeight = 9; % cm
FigFontSize = 8; % pt
type = 'c_lambda';
up_low_boundary = 'yes';
figure
MASWaves_plot_dispersion_curve(f_curve0,c_curve0,lambda_curve0,...
     f_curve0_up,c_curve0_up,lambda_curve0_up,f_curve0_low,c_curve0_low,...
     lambda_curve0_low,type,up_low_boundary,FigWidth,FigHeight,FigFontSize)
 

%% Inversion
c_test_min = 0; % m/s
c_test_max = 2500; % m/+s
delta_c_test = 50; % m/s
c_test = c_test_min:delta_c_test:c_test_max; % m/s

% Layer parameters
n = 10;  % number of layers used in the TLBO optimization
alpha = [4000 4000	4000	4000	4000	4000	4000 4000 4000 4000 4000]; % m/s  P-wave velocity
  h = [2.0	2.5	1	0.1	1.2	0.8	0.4	0.5	3.8	1.5	0.0]; % m soil layer thickness
beta = [150.0	150.0	1060.4	188.6	188.6	2000.0	1850.0	1442.3	310.5	150.0	1809.5]; % m/s Shear wave velocity
rho = [2000 2000 2000 2000 2000 2000 2000 2000 2000 2000 2000 2000 2000 2000 2000]; % kg/m^3 soil layer density
up_low_boundary = 'yes';
[c_t,lambda_t] = MASWaves_theoretical_dispersion_curve...
    (c_test,lambda_curve0,h,alpha,beta,rho,n);

up_low_boundary = 'yes';
FigWidth = 8; % cm
FigHeight = 8; % cm
FigFontSize = 12; % pt
figure
MASWaves_plot_theor_exp_dispersion_curves(c_t,lambda_t,...
    c_curve0,lambda_curve0,c_curve0_up,lambda_curve0_up,...
    c_curve0_low,lambda_curve0_low,up_low_boundary,...
    FigWidth,FigHeight,FigFontSize)

e = MASWaves_misfit(c_t,c_curve0);

f_curvet = f_curve0';

up_low_boundary = 'yes';
FigWidth = 16; % cm
FigHeight = 8; % cm
FigFontSize = 12; % pt
figure

MASWaves_plot_inversion_results_one_iteation(c_t,f_curvet,...
    c_curve0,f_curve0,c_curve0_up,f_curve0_up,c_curve0_low,...
    f_curve0_low,n,beta,h,e,up_low_boundary,FigWidth,FigHeight,FigFontSize)


save Result.mat
clock
toc

