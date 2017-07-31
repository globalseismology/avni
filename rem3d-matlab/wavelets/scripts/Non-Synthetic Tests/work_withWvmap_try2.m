%%Work with a Map of Wv Coeffs, test results from African topography. 
Choice_1 =0;
 ORIG = 'Direct'
setenv('IFILES','/home/moulik/Software/fjsimons-MMASS-v1.0.1/DATA')
colperc = [5 95];
Get_Edgy = 0;
N = 7; Jmax = 4; face = 3; J = 4;
defval('L',ceil(2^(N+1)))
defval('colmap','kelicol')
xbound =[40 90];
ybound = [60 110];
[ xpts,ypts ] = Get_Interior_PolygonPts_on_CubedSphere( N,[ybound(1) ybound(1) ybound(2) ybound(2)],[xbound(1) xbound(2) xbound(2) xbound(1)] );

[vwlev,vwlevs] = cube2scale(N,[Jmax Jmax],1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load Real Structure- In this case, Topography. 
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


African_Cubed_Chunk = v(:,:,face);
Reconstruct_Africa_Direct = angularD4WT(vw,[J J],[1 1],'inverse',1);
vw=angularD4WT(v,[Jmax Jmax],[1 1],'forward',1);

dax=prctile(African_Cubed_Chunk(:),colperc);
subplot(3,3,1)
h=imagefnan([1 1],[2^N 2^N],v(:,:,face),colmap,dax,[],[],0);
title('original anomaly')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

all_wavelets = angularD4WT(African_Cubed_Chunk,[Jmax Jmax],[1 1],'forward',1);
subplot(3,3,2)
h=imagefnan([1 1],[2^N 2^N],all_wavelets,colmap,dax,[],[],0);
title('wavelet coefficients of original anomaly')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MASK ORIGINAL REGION
if ORIG == 'Direct'
Structure_of_Interest = zeros(2^N,2^N,6);
    for i = 1:length(xpts)
        
        Structure_of_Interest(xpts(i),ypts(i),face) = Structure_of_Interest(xpts(i),ypts(i),face) + African_Cubed_Chunk(xpts(i),ypts(i));
    end
else
Structure_of_Interest = zeros(2^N,2^N,6);
    for i = 1:length(xpts)
        
        Structure_of_Interest(xpts(i),ypts(i),face) = Structure_of_Interest(xpts(i),ypts(i),face) + Reconstruct_Africa_Direct(xpts(i),ypts(i),face);
    end  
end
    subplot(3,3,3)
    aa=imagefnan([1 1],[2^N 2^N],Structure_of_Interest(:,:,face),colmap,dax,[],[],0);

    title('Masked, Original Geographical Region of Interest')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now,find locations of wavelet coefficients. 
if Choice_1 == 1
%%%%%%%%%%%%5
file = load('ForwardwaveletTry2.N7.J4.D4.mat');
Wv_Coeffs_Map = file.All_map;

Wavelet_Builder = zeros(128,128);
Zero_mat = zeros(128,128);
for i = 1:length(xpts)
    for j = 1:length(ypts)
        Temp_List = find(full(Wv_Coeffs_Map(xpts(i)).map{ypts(j)})~= 0);
        Zero_mat(Temp_List) = 1;
        Wavelet_Builder = Wavelet_Builder + Zero_mat;
        Zero_mat = zeros(128,128);
    end
    100*i/length(xpts)
end


else
%%%%%
vwt_inv = angularD4WT(Structure_of_Interest(:,:,face),[Jmax Jmax],[1 1],'forward',1);
for i = 1:length(vwt_inv)
    if vwt_inv(i) ~=0
        vwt_inv(i) = 1;
    end
end

Wavelet_Builder = vwt_inv;
end
subplot(3,3,4)
h=imagefnan([1 1],[2^N 2^N],Wavelet_Builder,colmap,dax,[],[],0);
title('Locations of wavelets formed by identifying wavelet maps corresponding to DV anomalies')



%%%%%%%%%%%%%%%%%%%%%%%%%

Orig_Wavelet_Transform = angularD4WT(African_Cubed_Chunk,[Jmax Jmax],[1 1],'forward',1);
Wavelet_locs = find(Wavelet_Builder ~= 0);
The_Big_One = zeros(2^N,2^N);
The_Big_One(Wavelet_locs) = Orig_Wavelet_Transform(Wavelet_locs);
subplot(3,3,5)
GiveItBackYouMonster = angularD4WT(The_Big_One,[Jmax Jmax],[1 1],'inverse',1);
hhhhh=imagefnan([1 1],[2^N 2^N],GiveItBackYouMonster,colmap,dax,[],[],0);

subplot(3,3,6)
resid = GiveItBackYouMonster - Structure_of_Interest(:,:,face);
contourf(resid)
hold on
plot([xbound(1) xbound(2) xbound(2) xbound(1) xbound(1)],[ybound(1) ybound(1) ybound(2) ybound(2) ybound(1)],'r')
colorbar
caxis([-5000 1])
