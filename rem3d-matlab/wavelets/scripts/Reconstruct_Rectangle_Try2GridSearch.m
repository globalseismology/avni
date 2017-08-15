%Try 2: Reconstruct Rectangle by finding all the wavelets which contribute
%anything to its location!
face = 3;
N = 7;
Jmax = 4;
Threshold = 0;
%%%Define Rectangle Locations

%Make the Rectangle!


Desired_Index_dim2 = [64 64 64 64 64 64 64 65 66 67 68 69 70 71 72 73 74 75 76 77 77 77 77 77 77 77 76 75 74 73 72 71 70 69 68 67 66 65];
Desired_Index_dim1 = [64 65 66 67 68 69 70 70 70 70 70 70 70 70 70 70 70 70 70 70 69 68 67 66 65 64 64 64 64 64 64 64 64 64 64 64 64 64];
Desired_Coefficient = ones(38,1);


[vwlev,vwlevs] = cube2scale(N,[Jmax Jmax],1);
TransformMe = zeros(2^N,2^N,6);
for i = 1:length(Desired_Index_dim1)
    %disp(i)
    TransformMe(Desired_Index_dim1(i),Desired_Index_dim2(i),face) = Desired_Coefficient(i);
end

figure
subplot(3,2,1)

i=imagefnan([1 1],[2^N 2^N],TransformMe(:,:,face));
title('synthetic velocity anomaly')

vws = angularD4WT(TransformMe,[Jmax Jmax],[1 1],'forward',1);
subplot(3,2,2)
j=imagefnan([1 1],[2^N 2^N],vws(:,:,face));
title('Wavelet Coefficients')

%%%%%%%%%%%%%%%%
Wavelet_Perturbation_Map = load('invwavelet.N.J.D4.mat');
Wavelet_Perturbation_Map = Wavelet_Perturbation_Map.All_map;
%Find Locations
Locations_of_Interest = find(TransformMe(:,:,face) ~= 999999);
%Preallocate Dv field with which we will reconstruct anomalies. 
Reconstruct_Dv_Field = zeros(2^N,2^N,6);
%


for i = 1:length(Locations_of_Interest)
    Curr_Location = Locations_of_Interest(i);
    for ii = 1:2^N
        for jj = 1:2^N
            if abs(Wavelet_Perturbation_Map(ii).map{jj}(Curr_Location)) > Threshold
                Reconstruct_Dv_Field(:,:,face)  = full(Wavelet_Perturbation_Map(ii).map{jj})*vws(ii,jj,face) + Reconstruct_Dv_Field(:,:,face);
            end
        end
    end
    (i/length(Locations_of_Interest))*100    %Progress bar in command line
end

subplot(3,2,3)
j=imagefnan([1 1],[2^N 2^N],Reconstruct_Dv_Field(:,:,face));
title(['Grid Search Reconstruction: Tol = ' num2str(Threshold)])
