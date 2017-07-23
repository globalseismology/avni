%Script to plot a given wavelet or wavelets on the sphere. 

%Essentially, can hard build a synthetic model element-by-element
%by scratch with this. Just be aware of where the scales are meant to be on
%the cubed sphere. 

%Anant Hariharan
face = 3;
Desired_Index_dim1 = [64 64 64 64 65 66 67 67 67 67 66 65];
Desired_Index_dim2 = [64 65 66 67 67 67 67 66 65 64 64 64];
Desired_Coefficient = [1 1 1 1 1 1 1 1 1 1 1 1];

N = 7;
Jmax = 4;

[vwlev,vwlevs] = cube2scale(N,[Jmax Jmax],1);
TransformMe = zeros(2^N,2^N,6);


for i = 1:length(Desired_Index_dim1)
    disp(i)
    TransformMe(Desired_Index_dim1(i),Desired_Index_dim2(i),face) = Desired_Coefficient(i);
end

reconstructed_synthDv = angularD4WT(TransformMe,[Jmax Jmax],[1 1],'inverse',1);

figure
subplot(1,2,1)

h=imagefnan([1 1],[2^N 2^N],reconstructed_synthDv(:,:,face));
title('reconstruction of 4 scale 4 wavelets')

subplot(1,2,2)

h=imagefnan([1 1],[2^N 2^N],TransformMe(:,:,face));
title('position of wavelets')




%%%%Trying to figure out how the hell this thing works. Can comment
%%%%everything below this line out if you want- probably won't be useful

%Anant Hariharan
%%% Cross
face = 3;
Desired_Index_dim1 = [63 64 65 63 65 62 66];
Desired_Index_dim2 = [65 64 65 63 63 62 66];
Desired_Coefficient = [100 100 100 100 100 100 100];

%%% Square

face = 3;
Desired_Index_dim1 = [64 64 64 64 65 66 67 67 67 67 66 65];
Desired_Index_dim2 = [64 65 66 67 67 67 67 66 65 64 64 64];
Desired_Coefficient = [1 1 1 1 1 1 1 1 1 1 1 1];
N = 7;
Jmax = 4;
[vwlev,vwlevs] = cube2scale(N,[Jmax Jmax],1);
TransformMe = zeros(2^N,2^N,6);
for i = 1:length(Desired_Index_dim1)
    disp(i)
    TransformMe(Desired_Index_dim1(i),Desired_Index_dim2(i),face) = Desired_Coefficient(i);
end
vws = angularD4WT(TransformMe,[Jmax Jmax],[1 1],'forward',1);

figure
subplot(1,2,1)

h=imagefnan([1 1],[2^N 2^N],vws(:,:,face));
title('wavelet coefficients describing anomaly')

subplot(1,2,2)

h=imagefnan([1 1],[2^N 2^N],TransformMe(:,:,face));
title('synthetic velocity anomaly')



