%%Test out the D6 wavelet inversion.
N = 9;
Jmax = 4;
Transform_Grid = zeros(2^N,2^N);
Transform_Grid(1) = 1;
%Transform_Grid(2) = 1;
recovered = angularD6WT(Transform_Grid,[Jmax Jmax],[1 1],'inverse',1);

figure
subplot(1,2,1)

h=imagefnan([1 1],[2^N 2^N],recovered);
title('reconstruction of a D6 wavelet')

subplot(1,2,2)

h=imagefnan([1 1],[2^N 2^N],Transform_Grid);
title('position of wavelets')