%% Geographic Thresholding with wavelets: MATLAB version.
%% Uses ONLY the .mat files set up in the rem3d/rem3d/files directory  
%% to extract a specific geographic region of a tomographic model in the wavelet basis. 

% Anant Hariharan
clear all
clc
close all
addpath(genpath('/home/anant/Software/rem3d/rem3d-matlab/wavelets'))

%% Load and Read Files... N  = 7, Superchunks. 
Grids = load('/home/anant/Software/rem3d/rem3d/files/SC_Grid_N7_Jmax4_EulerConfig1.mat');
WaveletMaps_DV_Basis = load('/home/anant/Software/rem3d/rem3d/files/SC_Mmaps.N7.J4.D4.mat');
WaveletCoeffstoMakeDV = load('/home/anant/Software/rem3d/rem3d/files/SC_WCoeffs.N7.J4.D4.mat');
ModelFile = load('SC_MIT_P08.N7.Jmax4.EulerConfig1.mat');

%% Set Params
% Note, Africa is on face 3. 
N=7;
LonMin = 5;
LonMax = 40;
LatMin = -14;
LatMax = 15;
Desired_Scales = [5]; %Not using this right now, but I should be able to do so easily. 

re = 6371;
depth=100;
desired_radius = re - depth;

%% Do Thresholding using the "Grids" file as a geographic reference. This example only does rectangular.
Geog_Indices = find(Grids.lon > LonMin & Grids.lon < LonMax & Grids.lat > LatMin & Grids.lat<LatMax);

%% Now, find the model coefficients that correspond to the points you desire. Put in a structure. 
Data.Indices = Geog_Indices;

% First, extract the coefficients for the right depth. 
ModelwvAtDepth = ModelFile.wvcoeffs(find(ModelFile.depth == desired_radius));
ModelveloAtDepth = ModelFile.model(find(ModelFile.depth == desired_radius));
%Data.velocoeffs

for i = 1:length(Data.Indices)
    Data.Wvcoeffs(i) = ModelwvAtDepth(Data.Indices(i));
    Data.Velocoeffs(i) = ModelveloAtDepth(Data.Indices(i));
end

%Edit the coeffs to be on the right face- can do this with rem better,
%probably. but hardwiring for this example for simplicity.

Data.Indices = Data.Indices-2*(2^(2*N));

%% Now, make a map of the wavelet coefficients that correspond to the thresholded region. 
%% Note: This is essentially a 'selective' forward transform into the wavelet domain. 

Reconstructed_Wavelet_Map = zeros(2^N,2^N);

for i = 1:length(Data.Indices)
    CurrMap = Data.Velocoeffs(i)*full(WaveletCoeffstoMakeDV.Me(Data.Indices(i)).dv);
    Reconstructed_Wavelet_Map = Reconstructed_Wavelet_Map+CurrMap;
end

figure
subplot(4,1,1)
h=imagefnan([1 1],[2^N 2^N],Reconstructed_Wavelet_Map)
title('Map of the wavelets used in the reconstruction')

%% Now do the inverse transform! 
%%

Reconstructed_VeloMap = zeros(2^N,2^N);
for i = 1:2^(2*N)
    CurrMap = Reconstructed_Wavelet_Map(i)*full(WaveletMaps_DV_Basis.Me(i).dv);
    Reconstructed_VeloMap = Reconstructed_VeloMap+CurrMap;
end

subplot(4,1,2)
contourf(Reconstructed_VeloMap)
hold on
[I,J] = ind2sub(size(Reconstructed_VeloMap),Geog_Indices)
scatter(I,J)
title('Thresholded Region')

%% Compare to the entire sixth face. 
Reconstructed_Model = zeros(2^N,2^N);
for i = 1:2^(2*N)
    CurrMap = ModelwvAtDepth(2*98304/6+i)*full(WaveletMaps_DV_Basis.Me(i).dv);
    Reconstructed_Model = Reconstructed_Model+CurrMap;
end
subplot(4,1,3)
contourf(Reconstructed_Model)
title('Whole face')

title('Full Face on Sphere')
% %% Plot on the face of a cartesian grid. 
% figure
% %Iterate through every value in the grid. Find the lat/lon. 
% for i = 1:2^(2*N)
%     scatter(Grids.lon(2*98304/6+i),Grids.lat(2*98304/6+i),5,Reconstructed_VeloMap(i))
%     hold on
% end
% grid on
% box on
% 

% Residual
Residual = Reconstructed_Model-Reconstructed_VeloMap;
subplot(4,1,4)
contourf(Residual)
title('Residual with Face and Thresholded Region')

