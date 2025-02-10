clear all
close all
clc
tic
%% DATA INPUT 
% Filename = 'shot_5mo_2dx_2.dat';
Filename = 'Shot4_4_0.3s.dat';
% Filename = 'Shot4_0.3s.dat';
HeaderLines = 7;
fs = 4000; % Hz
N = 6;
x1 = 5; % m
dx = 1; % m
Direction = 'forward';

[u,t,Tmax,L,x] = MASWaves_read_data(Filename,HeaderLines,fs,N,dx,x1,Direction);

%% PLOT DATA 
du = 1/75;
FigWidth = 6; % cm
FigHeight = 8; % cm
FigFontSize = 8; % pt

figure
MASWaves_plot_data(u,N,dx,x1,L,t,Tmax,du,FigWidth,FigHeight,FigFontSize)

%% DISPERSION 
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
% select = 'numbers';
select = 'both';
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
n = 10;
% alpha = [600 600	600	600	600	600	600 600 600	600	600	600	600	600]; % m/s  AL
alpha = [4000 4000	4000	4000	4000	4000	4000 4000 4000 4000 4000]; % m/s  shot1
% alpha = [5000 5000	5000	5000	5000	5000	5000 5000 5000 5000 5000 5000 5000 5000 5000]; % m/s
%   h = [2	2	3	1	1	3	7]; % model3
% h = [3.2	3.8	0.5	1.3	2.3	4.0	10.0]; % shot1  
%   h = [1	3.4	0.5	0.4	2.7	3.3	1.6	1.3	3.1	0.2]; % m shot4_3
  h = [2.0	2.65	1	0.1	2.5	0.8	2.4	0.5	3.8	1.5	0.0]; % m shot4_3
%  h = [ 4.3	1.0	1.0	1.2	1.24	0.39	1.29	0.31	0.69	0.76	0.71	1.84	0.48	1.69] 
%  h = [1.01	1.67	0.41	0.53	0.71	0.71	0.78	0.46	1.29	1.74	0.62	0.50	1.20	1.17]; % m shot4_4
% h = [3.4	2.1	0.3	3.3	4.0	4.0	3.0]; % m shot6
%  h = [0.58	0.76	0.39	0.29	0.46	0.68	1.04	0.90	0.18]; % m shot7
% beta = [154	1010	120	111	2297	1071	745	1656	163	1104	2609	891	2444	2810]; % m/s AL1
% beta = [166	135	120	278	299	300	122]; % m/s AL1
% beta = [164.25	555.41	163.63	116.70	799.48	132.04	509.78	106.24	100.00]; % m/s shot7
%  beta = [200	200	682	789	200	200	800]; % m/s shot6
% beta = [200	228	200	200	200	214	796]; % m/s shot5
% beta = [200	201	200	786	800	800	800]; % m/s shot1
beta = [150.0	150.0	1660.4	188.6	188.6	200.0	200.0	1442.3	310.5	150.0	1809.5]; % m/s Shot4_4
% beta = [213.35	101.05	176.40	246.93	292.20	278.48	222.79	163.52	154.82	136.34	127.93	281.84	195.76	252.46]; % m/s shot4_3
% beta = [163	131	154	300	298	299	120]; % m/s AL4
% rho = [2000 2000 2000 2000 2000 2000 2000 2000 2000 2000 2000]; % kg/m^3
rho = [2000 2000 2000 2000 2000 2000 2000 2000 2000 2000 2000 2000 2000 2000 2000]; % kg/m^3
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
% MASWaves_plot_inversion_results_one_iteation(c_t,lambda_t,...
%     c_curve0,lambda_curve0,c_curve0_up,lambda_curve0_up,c_curve0_low,...
%     lambda_curve0_low,n,beta,h,e,up_low_boundary,FigWidth,FigHeight,FigFontSize)
MASWaves_plot_inversion_results_one_iteation(c_t,f_curvet,...
    c_curve0,f_curve0,c_curve0_up,f_curve0_up,c_curve0_low,...
    f_curve0_low,n,beta,h,e,up_low_boundary,FigWidth,FigHeight,FigFontSize)


save Result.mat
clock
toc

