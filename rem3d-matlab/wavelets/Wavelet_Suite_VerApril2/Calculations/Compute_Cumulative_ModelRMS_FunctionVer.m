function [Cumulative_RMS_List] = Compute_Cumulative_ModelRMS_FunctionVer(Model1,Geographic_Lon,Geographic_Lat,N,Jmax)
% This function runs through the different expansions of a model in the D4
% wavelet domain, calculating RMS at each one by iteratively adding in new wavelet
% coefficients that belong to each scale, and returns the list as a
% vector of RMSs. 
% Anant Hariharan, 3/3/18

% Standard assumptions, including model on cubed sphere already zeroed out, etc. 
Cumulative_RMS_List = zeros(1,6);
Scale_List = Jmax+1:-1:1;
Grids = load('/home/anant/Software/rem3d/rem3d/files/Csph_Grid_N8_Jmax5_EulerConfig1.mat');
% Locations of wavelet, i.e. scaling functions and such. 
vwlev1=cube2scale(N,[Jmax Jmax]+1,1,1);

% First threshold over the geographic region. 
[Model1_WaveletOut,Model1_PixelOut] = Threshold_Model_Geographic(N,Grids,Geographic_Lon,Geographic_Lat,Model1,Jmax);

for llll = 1:length(Scale_List)
llll
Desired_Scales = Scale_List(1:llll);
Reconstructed_ModelWavelet_Map_ScaleInspected = zeros(2^N,2^N);

for ii = 1:length(Desired_Scales)
    curr_scale = Desired_Scales(ii);
Scale_check = find(vwlev1 == curr_scale);
Reconstructed_ModelWavelet_Map_ScaleInspected(Scale_check) = Model1_WaveletOut(Scale_check);

end

%% Now do the inverse transform! 
%%
Reconstructed_ModelVeloMap = angularD4WT(Reconstructed_ModelWavelet_Map_ScaleInspected,[Jmax Jmax],[1 1],'inverse',1);

Curr_ModelRMS = sqrt(mean(mean(Reconstructed_ModelVeloMap.^2)));

Cumulative_RMS_List(llll)= Curr_ModelRMS;

end
end




















