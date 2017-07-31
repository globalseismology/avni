%%%Make Comparison Subplots of Wavelet Results

results = load('loris5_GN_0023_7_4_D4_1_1.mat');
vw = results.vw;
figure(1)
plotoncube(vw,'2D')
vwrec = angularD4WT(vw,[4 4],[1 1],'inverse',1);
figure(2)
plotoncube(vwrec,'2D')

subplot(1,3,1); 
make2d_map_plot(vwrec);
title('Inversion and Reconstruction at all Scales');
xlabel('longitude')
ylabel('latitude')

%reconstruct for all scale 3 structure.
scale4indices = getkeepindex(vw,7,4,4);
scale4wavelets = zero_wavelets(vw,scale4indices);

vwrecscale4 = angularD4WT(scale4wavelets,[4 4],[1 1],'inverse',1);
figure(3)
plotoncube(vwrecscale4,'2D')
subplot(1,3,2); 
make2d_map_plot(vwrecscale4);

title('Inversion and Reconstruction at Scale 4');
xlabel('longitude')
ylabel('latitude')

%reconstruct for all scale 4 structure and NA structure
scale4indicesNA = getkeepindex(vw,7,4,4,-90,90,-180,0);
%latmin,latmax,lomin,lomax getkeepindex(vw,N,Jmin,Jmax,latmin,latmax,lomin,lomax
scale4waveletsNA = selected_wavelets2dv(vw,scale4indicesNA);
vwrecscale4NA = angularD4WT(scale4waveletsNA,[4 4],[1 1],'inverse',1);

subplot(1,3,3); 
make2d_map_plot(vwrecscale4NA);
caxis([-3 3])
title('Inversion and Reconstruction at Scale 4 for wavelets on The Americas');
xlabel('longitude')
ylabel('latitude')