% Set the defaults
defval('wav','D4')
defval('spp',1)
defval('N',7)
defval('J',4)
defval('L',ceil(2^(N+1)))
defval('precon',[1 1]*0)
defval('iface',3)
defval('colmap','kelicol')
setenv('IFILES','/home/moulik/Software/fjsimons-MMASS-v1.0.1/DATA')
% Color map saturation as percentiles
colperc=[5 95];
% For SPIE change it
colperc=[10 90];

% Truncation level as percentiles
defval('tperc',85);
% More percentiles for the tick marks on the color bar
ticperc=[50];

% Load or make the data, that is, the topography
fname=fullfile(getenv('IFILES'),'EARTHMODELS','CONSTANTS',...
	       sprintf('loris2_%i_%i.mat',N,L));

if exist(fname,'file')==2
  load(fname)
else
  % Load Earth's topography
  lmcosi=rindeks(fralmanac('GTM3AR','SHM'),1:addmup(L));
  
  % Perform the transform with standard inputs
  v=plm2cube(lmcosi,N);

  % Save the data for the next time you run this
  save(fname,'v','N','L')
end
colperc = [5 95];

vw=angularD4WT(v,[J J],[1 1],'forward',1);

indices_minusscale2 = getkeepindex(vw,7,4,[1 3 4])

vw_minusscale2 = zero_wavelets(vw,indices_minusscale2);

reconstructed_topo = angularD4WT(vw_minusscale2,[J J],precon,'inverse',1);

dax=prctile(reconstructed_topo(:),colperc);
h=imagefnan([1 1],[2^N 2^N],v(:,:,iface),colmap,dax,[],[],0);
hold on
% Plot the wavelet grid
  f=fridplotw(N,J);