%%%%
%TESTS 7-24

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

subplot(2,2,1)
j=imagefnan([1 1],[2^N 2^N],TransformMe(:,:,face));
title(['ORIGINAL'])

subplot(2,2,2)
h=imagefnan([1 1],[2^N 2^N],vws(:,:,face));
title(['wavelet coeffs of orig. struc'])

reconstruct = angularD4WT(vws,[Jmax Jmax],[1 1],'inverse',1);

subplot(2,2,3)
h=imagefnan([1 1],[2^N 2^N],reconstruct(:,:,face));
title(['Reconstruction from wavelet coeffs'])
