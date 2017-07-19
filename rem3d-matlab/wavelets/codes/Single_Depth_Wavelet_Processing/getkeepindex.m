function [ indicesout ] = getkeepindex(vw,N,Jup,J_Interest,latmin,latmax,lomin,lomax,thresholdpercentile,eo)
% Given a desired method of thresholding, we separate out the wavelet
% indices on the cubed sphere that correspond to a desired range- this range may be 
% spatially constrained or constrained by a desired scale range, or percentile 

% Significant information from Frederik Simons' Loris5 codes. 

% - Anant Hariharan

defval('N',7);
defval('Jup',4);
defval('Jmin',4);
defval('Jmax',4);
defval('latmin',-90);
defval('latmax',90);
defval('lomin',-180);
defval('lomax',180);
defval('thresholdpercentile',100);
defval('eo',0);


numvals = 6*2^(2*N);
%Get indices map that link points on the cubed sphere to a specific scale. 
[vwlev,vwlevs]=cube2scale(N,[Jup Jup],1);

%Preallocate using resolution parameter
scales_passthetest = zeros(2^N,2^N,6);
%%%% First, separate out by scale!
Scales_Interest = J_Interest;
for i = 1:length(Scales_Interest)
    jndex = Scales_Interest(i);
    itshere=[vwlevs==jndex];
    scales_passthetest = scales_passthetest+itshere;
end


%%%% Separate out by thesholding percentile.
perc_passthetest=vw<prctile(vw(:), thresholdpercentile);


%%%% Separate out lat/long!
[x,y,z]=cube2sphere(N,[],[],[],eo);
% Each run takes care of one chunk on the cubed sphere. 
latlon_passthetest = zeros(2^N,2^N,6);
counter = 0;
for in=1:6
   
   % Transform these to longitude and latitude
   [phi,piminth,r]=cart2sph(x(:,:,in),y(:,:,in),z(:,:,in));
   lon=phi(:)*180/pi; lat=piminth(:)*180/pi;
   for i = 1:length(lon)
       counter = counter + 1;
       %Does it match the lat/lon bounds?
       if lon(i) >= lomin && lon (i) <= lomax && lat(i) >= latmin && lat (i) <= latmax
           latlon_passthetest(counter) = 1;
       end
       
   end
   
end

%indicesout = (scales_passthetest  && latlon_passthetest)>0) == perc_passthetest == ones(2^N,2^N,6);

 indicesout = zeros(2^N,2^N,6);
 
 for i = 1:numvals
     if scales_passthetest(i)>0 && latlon_passthetest(i)>0 && perc_passthetest(i)>0
         indicesout(i)=1;
     end
 end



end

