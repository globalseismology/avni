function [ output_args ] = make2d_map_plot( vw,N,J,eo,size )
% Makes map plot of cubed sphere data.
defval('N',7);
defval('J',4);
defval('eo',0);
defval('size',10);
numvals = 6*2^(2*N);
lfin =N;
[x,y,z]=cube2sphere(lfin,[],[],[],eo);
megalon = [];
megalat = [];

for in=1:6
   
   % Transform these to longitude and latitude
   [phi,piminth,r]=cart2sph(x(:,:,in),y(:,:,in),z(:,:,in));
   lon=phi(:)*180/pi; lat=piminth(:)*180/pi;
   megalon = [megalon; lon]; 
   megalat = [megalat; lat];
   
end

coeffs = [];
for i = 1:numvals
    coeffs = [coeffs vw(i)];
end

scatter(megalon,megalat,size,coeffs,'filled')
colorbar
%caxis([-0.05 0.05])
%title(['Absolute Value of Wavelet Coefficients for Depth =' dep 'km'])
xlabel('longitude')
ylabel('latitude')
end

