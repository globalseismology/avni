function [ maxdeg,WavelengthKm] = Scale2Wavelength_Depth( N,Scale,Jmax,depkm )
% Calculate the dominant wavelength and degree of an equivalent wavelet by interpolating on
% grid, converting to SH deg, then using the Jean Relation
% (see treatise on geophysics, page 155) to get dominant wavelength.

% Do this given the resolution parameter,and a given depth(km) within the Earth.  
% For daubechies wavelet.

% Anant Hariharan, but lines 19-24 is from loris6, etc.
Re=6371;
Transform_Me = zeros(2^N,2^N,6);
wavelet_locs = get_wavelet_pos( 2^(N-1),2^(N-1),N,Scale,1,0);
[vwlev,vwlevs] = cube2scale(N,[Jmax Jmax],1);
pts_interest = find(vwlevs == Scale);
Transform_Me(pts_interest(1)) = 1;
vws = angularD4WT(Transform_Me,[Jmax Jmax],[1 1],'inverse',1);


[lons,lats,vd,lon,lat]=cube2lola(vws); 
vd(isnan(vd))=0;
% Figure out the dominant deg.
[sdl,l]=plm2spec(xyz2plm_hardfix(vd));
%semilogy(l,sdl)
[a,maxdeg]=max(sdl);

%%%Now, use the Jean Relation to get actual cartesian wavelength from degree.

WavelengthKm = 2*pi*(Re-depkm)/(sqrt(maxdeg*(maxdeg+1)));

end