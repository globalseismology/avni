%Script to plot a given wavelet or wavelets on the sphere. 

%Essentially, can hard build a synthetic model element-by-element
%by scratch with this. Just be aware of where the scales are meant to be on
%the cubed sphere. 

%Anant Hariharan
face = 3;
Desired_Index_dim1 = [ 1 16 ];
Desired_Index_dim2 = [ 1 1];
Desired_Coefficient =[100 100];

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


figure
subplot(1,2,1)
title('reconstruction of 4 scale 4 wavelets')
h=imagefnan([1 1],[2^N 2^N],reconstructed_synthDv(:,:,face));subplot(1,2,2)
subplot(1,2,2)
title('position of wavelets')
h=imagefnan([1 1],[2^N 2^N],TransformMe(:,:,face));