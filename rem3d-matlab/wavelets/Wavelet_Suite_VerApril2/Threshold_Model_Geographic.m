function [WaveletOut,PixelOut] = Threshold_Model_Geographic(N,Grids,PolyLon,PolyLat,Csph_Face,Jmax)
%% Illustrate the Workflow of Geographic Thresholding
[in,on] = inpolygon( Grids.lon,Grids.lat,PolyLon,PolyLat);   % Logical Matrix
inon = in | on;     % Combine ?in? And ?on?
idx = find(inon(:));  
Geog_Indices = idx;

%%%%%%%%%%%%%%%%%%%%%% Step 1: Mask Structure of Interest
Masked_Model = zeros(2^N,2^N);
Masked_Model(idx) = Csph_Face(idx);
 
%%%%%%%%%%%%%%%%%%%%%% Step 2: Get the locations of Wavelet Coefficients of Structure of
%%%%%%%%%%%%%%%%%%%%%% Interest
Model_Mask_vwt_inv = angularD4WT(Masked_Model(:,:),[Jmax Jmax],[1 1],'forward',1);
Model_Wavelet_locs = find(Model_Mask_vwt_inv ~= 0);

%%%%%%%%%%%%%%%%%%%%%% Step 3: Get the Original Wavelet Coefficients at
%%%%%%%%%%%%%%%%%%%%%% these locations
 
Orig_Wavelet_Transform = angularD4WT(Csph_Face,[Jmax Jmax],[1 1],'forward',1);
Final_Wavelet_Coeff_Map = zeros(2^N,2^N);
Final_Wavelet_Coeff_Map(Model_Wavelet_locs) = Orig_Wavelet_Transform(Model_Wavelet_locs);
 
%%%%%%%%%%%%%%%%%%%%%% Step 4: Do the inverse transform at these locations.
Thresholded_Out = angularD4WT(Final_Wavelet_Coeff_Map,[Jmax Jmax],[1 1],'inverse',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PixelOut = Thresholded_Out;
WaveletOut = Final_Wavelet_Coeff_Map;
end

