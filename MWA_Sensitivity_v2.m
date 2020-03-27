clear all;
close all;

%!============== VARIABLES ================!
%source pointing angles in degrees
phi0 = 0; %'phi' location of the source
theta0 = 0; %'theta' location of the source
gridpoint = 0; %gridpoint number from sweetpoints
dPhi = 5; %phi resolution (degrees per angular step)
dTheta = 5; %theta resolution (degrees per angular step)
eff_rad = 0.95; %radiation efficiency

%!========= END OF VARIABLES =============!

T0 = 290; % Kelvin, ambient temperature
c0 = 299792458; %Speed of Light in Vacuum, m/s
kB = 1.380649e-23; %J.K-1 (Boltzmann constant)
M = 16; %number of elements

%=== noise parameters provided from Daniel ===
load('Data\EDA_LNA_DISO_170705'); 

freq_noise = noise_par_DISO.freq; % MHz
Tmin = noise_par_DISO.Tmin;
Gamma_opt = noise_par_DISO.gamma_opt_comp;
N = noise_par_DISO.N;

clear noise_par_DISO

%===calculate differential S-parameters ===
Spara = sparameters('Data\WORKIB LNA.S3P');
Spara = newref(Spara,100);

S = Spara.Parameters;
freq_DUT = Spara.Frequencies/1e6;

Sdd11(:,1) = 0.5 * (S(1,1,:) - S(2,1,:) - S(1,2,:) + S(2,2,:));
Sds12(:,1) = 1.0 / sqrt(2) * (S(1,3,:) - S(2,3,:));
Ssd21(:,1) = 1.0 / sqrt(2) * (S(3,1,:) - S(3,2,:));
Sss22(:,1) = S(3,3,:);

clearvars Spara S 

% %---firstly, perform interpolation to match data points---
Sdd11 = spline(freq_DUT,Sdd11,freq_noise);
Sds12 = spline(freq_DUT,Sds12,freq_noise);
Ssd21 = spline(freq_DUT,Ssd21,freq_noise);
Sss22 = spline(freq_DUT,Sss22,freq_noise);

%---now calculate the travelling noise wave correlation parameters 
c11 = kB * T0 *( 4 * N .* abs(1 - Sdd11 .* Gamma_opt).^2 ./ ...
        (1 - abs(Gamma_opt).^2) - Tmin / T0 .* (1 - abs(Sdd11).^2 ));
    
c22 = kB * T0 * abs(Ssd21).^2 .* (Tmin / T0 + ...
        4 * N .* abs(Gamma_opt).^2 ./ (1 - abs(Gamma_opt).^2));
    
c12 = -4 * kB * T0 * N .* (conj(Ssd21) .* conj(Gamma_opt)) ./ ...
      (1 - abs(Gamma_opt).^2) + Sdd11 ./ Ssd21 .* c22;  

%=== extract multiport scattering parameters from simulated data ===
S_ant_data = sparameters('Data\MWA_ARRAY_SParameter1.s16p'); 
freq_ant = S_ant_data.Frequencies/1e6;
freq_num = size(freq_ant,1);
S_ant = S_ant_data.Parameters;

clear S_ant_data

%=== interpolate to match the number of simulation frequency points ===
c11 = spline(freq_noise,c11,freq_ant);
c22 = spline(freq_noise,c22,freq_ant);
c12 = spline(freq_noise,c12,freq_ant);
 
Sdd11 = spline(freq_noise,Sdd11,freq_ant);
Sds12 = spline(freq_noise,Sds12,freq_ant);
Ssd21 = spline(freq_noise,Ssd21,freq_ant);
Sss22 = spline(freq_noise,Sss22,freq_ant);

Gamma_opt = spline(freq_noise,Gamma_opt,freq_ant);
N = spline(freq_noise,N,freq_ant);
Tmin = spline(freq_noise,Tmin,freq_ant);

%=== Read beamformer coefficients ===
csv_data = readtable('Data\MWA_sweet_spot_gridpoints.csv','HeaderLines',6);

%create complex beamformer coefficients, w
delays = table2array(csv_data(gridpoint+1, 4:end));

clear csv_data

amps = 1;
phases = 2 * pi * (freq_ant * 1e6) * (-delays) * 435e-12;
w = amps * exp(1i * phases);

Geff = zeros(freq_num, M);
Xi = zeros(freq_num, M);
Gamma_s = zeros(freq_num, M);
Gamma_eff = zeros(freq_num, M);

%=== Loop through all LNAs ===
for m = 1:M
    
    w_m = w(:,m);
    
    %=== Form interpolated S_lna matrix ===
    S_lna = zeros(M, M, freq_num);
    
    for j = 1:M
        if(j ~= m)
            S_lna(j, j, :) = Sdd11(:);
        end
    end
    
    %create Gamma_s
    GT = zeros(M, 1, freq_num);
    IM = eye(M);
    for f = 1:freq_num
        GT(:, :, f) = inv(IM(:, :) - S_ant(:, :, f) ...
                                   * S_lna(:, :, f)) * S_ant(:, m, f);
    end
    Gamma_s_m = squeeze(GT(m,:,:));
    
    %create bm
    c1 = zeros(M, 1);
    c1(m) = 1;
    
    bm = zeros(freq_num, M);
    for f = 1:freq_num
       bm(f, :) = inv(IM(:, :) - S_ant(:, :, f) * ...
                                 Sdd11(f)) * S_ant(:, :, f) * c1(:, 1);
    end
    
    %create kappa
    kappa = zeros(freq_num, 1);
    for i = 1:M
       if(i ~= m)
           for f = 1:freq_num
                kappa(f) = kappa(f) + ...
                    (w_m(f) * bm(f, i)) / (w_m(f) * bm(f, m)); 
           end 
       end
    end
    
    %create zeta
    zeta = (1 + kappa) ./ (1 + Gamma_s_m .* Sdd11 .* kappa);
    
    %create Gamma_eff
    Gamma_eff_m = zeta .* Gamma_s_m;
    
    %create Xi
    Xi_m = (1 - abs(Gamma_eff_m) .^ 2) ./ (1 - abs(Gamma_s_m) .^ 2);
    
    %create Geff
    Geff_m = (1 - abs(Gamma_s_m) .^ 2) ./ ...
              abs(1 - Gamma_s_m .* Sdd11) .^ 2 .* ...
              abs(Ssd21 .* w_m) .^ 2 .* ...
              abs(1 + Gamma_s_m .* Sdd11 .* kappa) .^ 2;
    
    Gamma_s(:, m) = Gamma_s_m(:);
    Gamma_eff(:, m) = Gamma_eff_m(:);
    Xi(:, m) = Xi_m(:);
    Geff(:, m) = Geff_m(:);
    
end

%calculate Gamma_opt_arr
Gamma_opt_arr = sum(Gamma_opt .* Gamma_s, 2) / M;

%calculate effective temperature
Teff = zeros(freq_num, M);

for m = 1:M
    Teff(:, m) = Tmin .* Xi(:, m) + 4 * T0 * N .* Xi(:, m) .* ...
                 abs(Gamma_eff(:, m) - Gamma_opt_arr).^2 ./ ...
                 ((1 - abs(Gamma_opt_arr) .^ 2) .* ...
                  (1 - abs(Gamma_eff(:, m)) .^ 2));
end

%calculate receiver noise temperature and average gain
Tout = sum(Teff .* Geff, 2);
Gavg = sum(Xi .* Geff, 2);
Trx = Tout ./ Gavg;

%calculate homogeneous sky temperature
lambda = c0 ./ (freq_ant * 1e6);

% Get Directivity Data 
%Read directivity saved in separate .csv files
phi_num = 73; %number of phi points
theta_num = 19; %number of theta points

MWA_Directivity = zeros(phi_num, theta_num, freq_num);

D_PATH = 'Data\Directivity\';

freq_dir = 50; %MHz (starting point)
for b=1:size(MWA_Directivity, 3)
  MWA_Directivity(:,:,b) = csvread(strcat(D_PATH, 'Directivity_', ...
                           num2str(freq_dir, '%.2f'), '.csv'), 0, 0);
  freq_dir = freq_dir + 1.28; %1.28 MHz step 
end

theta_array = [0:dTheta:90]';
phi_array = [0:dPhi:360]';

%find the closest position to phi0 location of the source
phi_pointed = 0;
phi_diff1 = abs(phi0 - phi_pointed);
phi_pos = 1;
for i = 2:phi_num
    phi_diff2 = abs(phi0 - phi_array(i));
    if(phi_diff2 < phi_diff1)
        phi_pointed = phi_array(i);
        phi_diff1 = phi_diff2;
        phi_pos = i;
    else
        break
    end
end

%find the closest position to theta0 location of the source
theta_pointed = 0;
theta_diff1 = abs(theta0 - theta_pointed);
theta_pos = 1;
for i = 2:theta_num
    theta_diff2 = abs(theta0 - theta_array(i));
    if(theta_diff2 < theta_diff1)
        theta_pointed = theta_array(i);
        theta_diff1 = theta_diff2;
        theta_pos = i;
    else
        break
    end
end

fprintf('The closest position to (phi,theta) = (%.2f,%.2f) is (%.2f,%.2f)\n', ...
        phi0, theta0, phi_pointed, theta_pointed);
    
%convert from dB to Linear
D0 = 10 .^ (squeeze(MWA_Directivity(phi_pos, theta_pos, :)) / 10); 
Aeff = lambda .^ 2 / (4 * pi) .* D0;

%=== calculate antenna noise temperature due to homogeneous sky ===
Tsky = 60 * lambda .^ 2.55;

%=== calculate Tsys, assuming radiation efficiency is 95% ===
Tsys = eff_rad * Tsky + (1 - eff_rad) * T0 + Trx;

%=== calculate Sensitivity ===
Sensitivity = Aeff ./ Tsys;

%=== calculate SEFD ====
SEFD = 2 * kB ./ Sensitivity * 1e26;

%=== export data for x and y polarizations of SEFD for comparison
load('Data\X_pol.mat');
load('Data\Y_pol.mat');

SEFD_Xpol(SEFD_Xpol == 1) = NaN;
SEFD_Xpol(SEFD_Xpol < 1e3) = NaN;
SEFD_Xpol(SEFD_Xpol > 1e7) = NaN;

SEFD_Ypol(SEFD_Ypol == 1) = NaN;
SEFD_Ypol(SEFD_Ypol < 1e3) = NaN;
SEFD_Ypol(SEFD_Ypol > 1e7) = NaN;

SEFD_Xpol = SEFD_Xpol(:);
SEFD_Ypol = SEFD_Ypol(:);
Freq_Xpol = Freq_Xpol(:);
Freq_Ypol = Freq_Ypol(:);

semilogy(Freq_Xpol,SEFD_Xpol,'.',...
         Freq_Ypol,SEFD_Ypol,'.',...
         freq_ant,SEFD);
legend("SEFD x-polarization", "SEFD y-polarization", "SEFD Calculated 0.9m");
xlabel("Frequency, MHz");
ylabel("System Equivalent Flux Density, Jansky");
title(strcat('Simulated and Measured SEFD Comparison at \phi=',...
      num2str(phi_pointed),' ^{o}, \theta = ', ...
      num2str(theta_pointed),'^{o} Pointing'))
xlim([50 350]);
grid on
