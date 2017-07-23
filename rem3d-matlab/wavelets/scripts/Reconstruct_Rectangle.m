%Reconstruction Tests!
%%%How accurately can we reconstruct a rectangle?

%Make the Rectangle!

face = 3;
Desired_Index_dim2 = [64 64 64 64 64 64 64 65 66 67 68 69 70 71 72 73 74 75 76 77 77 77 77 77 77 77 76 75 74 73 72 71 70 69 68 67 66 65];
Desired_Index_dim1 = [64 65 66 67 68 69 70 70 70 70 70 70 70 70 70 70 70 70 70 70 69 68 67 66 65 64 64 64 64 64 64 64 64 64 64 64 64 64];
Desired_Coefficient = ones(38,1);

N = 7;
Jmax = 4;
[vwlev,vwlevs] = cube2scale(N,[Jmax Jmax],1);
TransformMe = zeros(2^N,2^N,6);
for i = 1:length(Desired_Index_dim1)
    disp(i)
    TransformMe(Desired_Index_dim1(i),Desired_Index_dim2(i),face) = Desired_Coefficient(i);
end

vws = angularD4WT(TransformMe,[Jmax Jmax],[1 1],'forward',1);
% 
% figure
% subplot(1,2,1)
% 
% h=imagefnan([1 1],[2^N 2^N],vws(:,:,face));
% title('wavelet coefficients describing anomaly')
% 
% subplot(1,2,2)
% 
% h=imagefnan([1 1],[2^N 2^N],TransformMe(:,:,face));
% title('synthetic velocity anomaly')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Make half a rectangle
Indices_Less_Thanx71 = find(Desired_Index_dim2<71)
Desired_Index_dim2_halved = Desired_Index_dim2(Indices_Less_Thanx71)
Desired_Index_dim1_halved = Desired_Index_dim1(Indices_Less_Thanx71)
TransformMe2 = zeros(2^N,2^N,6);

for i = 1:length(Desired_Index_dim1_halved)
    disp(i)
    TransformMe2(Desired_Index_dim1_halved(i),Desired_Index_dim2_halved(i),face) = Desired_Coefficient(i);
end

%Detect Wavelets where we want them to be. 
[ wavelet_locs ] = get_wavelet_pos( Desired_Index_dim1_halved,Desired_Index_dim2_halved,N,Jmax,face,0)

%Now zero out wavelets at those locations. 
[ half_rectangle_wavelets ] = zero_wavelets(vws,wavelet_locs,'db4',7,1,4 )

%Reconstruct model using those wavelets
[half_rectangle_maybe] = angularD4WT(half_rectangle_wavelets,[Jmax Jmax],[1 1],'inverse',1);

figure
subplot(3,2,1)

h=imagefnan([1 1],[2^N 2^N],vws(:,:,face));
title('wavelet coefficients describing anomaly')

subplot(3,2,2)

i=imagefnan([1 1],[2^N 2^N],TransformMe(:,:,face));
title('synthetic velocity anomaly')

subplot(3,2,3)

j=imagefnan([1 1],[2^N 2^N],wavelet_locs(:,:,face));
title('Location of Wavelet Coefficients correponding to half rectangle')

subplot(3,2,4)

k=imagefnan([1 1],[2^N 2^N],TransformMe2(:,:,face));
title('Half Rectangle Anomaly')

%

subplot(3,2,5)

l=imagefnan([1 1],[2^N 2^N],half_rectangle_wavelets(:,:,face));
title('Zeroed Out Wavelets')

subplot(3,2,6)

m=imagefnan([1 1],[2^N 2^N],half_rectangle_maybe(:,:,face));
title('Reconstructed half rectangle? Anomaly')