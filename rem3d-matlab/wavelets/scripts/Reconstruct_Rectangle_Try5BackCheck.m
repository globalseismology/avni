
%Make the Rectangle!
%[ xpts,ypts ] = Get_Interior_PolygonPts_on_CubedSphere( N,xbounds,ybounds )
face = 3;
Desired_Index_dim2 = [64 64 64 64 64 64 64 65 66 67 68 69 70 71 72 73 74 75 76 77 77 77 77 77 77 77 76 75 74 73 72 71 70 69 68 67 66 65];
Desired_Index_dim1 = [64 65 66 67 68 69 70 70 70 70 70 70 70 70 70 70 70 70 70 70 69 68 67 66 65 64 64 64 64 64 64 64 64 64 64 64 64 64];
Desired_Coefficient = ones(2^25,1);

N = 8;
Jmax = 5;

[ Desired_Index_dim1,Desired_Index_dim2 ] = Get_Interior_PolygonPts_on_CubedSphere( N,[64 64 70 70],[64 70 70 64] );

[vwlev,vwlevs] = cube2scale(N,[Jmax Jmax],1);
TransformMe = zeros(2^N,2^N,6);
for i = 1:length(Desired_Index_dim1)
    %disp(i)
    TransformMe(Desired_Index_dim1(i),Desired_Index_dim2(i),face) = Desired_Coefficient(i);
end

vws = angularD4WT(TransformMe,[Jmax Jmax],[1 1],'forward',1);






[ xpts,ypts ] = Get_Interior_PolygonPts_on_CubedSphere( N,[64 64 70 70],[64 70 70 64] )

TransformMe2 = zeros(2^N,2^N,6);

for i = 1:length(xpts)
    disp(i)
    TransformMe2(xpts(i),ypts(i),face) = Desired_Coefficient(i);
end



vws_locs = angularD4WT(TransformMe2,[Jmax Jmax],[1 1],'inverse',1);


%Now keep wavelets at those locations. 
[ half_rectangle_wavelets ] = zero_wavelets(vws,vws_locs~=0,'db4',7,1,4 )

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

% j=imagefnan([1 1],[2^N 2^N],wavelet_locs(:,:,face));
% title('Location of Wavelet Coefficients correponding to half rectangle')

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