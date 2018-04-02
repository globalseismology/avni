%%Work with a Map of Wv Coeffs, test results from African topography. 

setenv('IFILES','/home/moulik/Software/fjsimons-MMASS-v1.0.1/DATA')
colperc = [5 95];
Get_Edgy = 0;
N = 7; Jmax = 4; face = 3; J = 4;
defval('L',ceil(2^(N+1)))
defval('colmap','kelicol')
xbound =[40 90];
ybound = [60 110];
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


dax=prctile(African_Cubed_Chunk(:),colperc);
subplot(3,3,1)
h=imagefnan([1 1],[2^N 2^N],v(:,:,face),colmap,dax,[],[],0);
title('original anomaly')

%Let's say the EARS is around x = 40:90, y = 60:120 on cubed sphere for N =
%7. Hard thresholding. 
[ xpts,ypts ] = Get_Interior_PolygonPts_on_CubedSphere( N,[ybound(1) ybound(1) ybound(2) ybound(2)],[xbound(1) xbound(2) xbound(2) xbound(1)] );
subplot(3,3,2)
plot(ypts,xpts,'ro')
title('points of interest')




Structure_of_Interest = zeros(2^N,2^N,6);
for i = 1:length(xpts)       
    Structure_of_Interest(xpts(i),ypts(i),face) = Structure_of_Interest(xpts(i),ypts(i),face) + African_Cubed_Chunk(xpts(i),ypts(i));
end

subplot(3,3,3)
aa=imagefnan([1 1],[2^N 2^N],Structure_of_Interest(:,:,face),colmap,dax,[],[],0);
title('Truncated, Original Geographical Region of Interest')





    

    
    
    
%%%%%%%%%%%%5
file = load('ForwardwaveletTry2.N7.J4.D4.mat')
Wv_Coeffs_Map = file.All_map;

Wavelet_Builder = zeros(128,128);
Zero_mat = zeros(128,128)
for i = 1:2^N
    for j = 1:2^N
%       Wavelet_Builder = Wavelet_Builder + Structure_of_Interest(i,j,face) * full(Wv_Coeffs_Map(i).map{j});
        Temp_List = find(full(Wv_Coeffs_Map(i).map{j})) ~= 0
        Zero_mat(Temp_List) = 1;
        Wavelet_Builder = Wavelet_Builder + Zero_mat;
        Zero_mat = zeros(128,128);
    end
end

subplot(3,3,4)
h=imagefnan([1 1],[2^N 2^N],Wavelet_Builder,colmap,dax,[],[],0);
title('wavelets formed by multiplying wavelet maps by DV anomalies')

Maybe_Res = angularD4WT(Wavelet_Builder,[Jmax Jmax],[1 1],'inverse',1);
subplot(3,3,5)
hhhhh=imagefnan([1 1],[2^N 2^N],Maybe_Res,colmap,dax,[],[],0);

%%%%%%
% Wvoneone_Map = full(Wv_Coeffs_Map(1).map{1});
% plotoncube(Wvoneone_Map,'2D')
% 
% 
% Orig_Wavelet_Transform = angularD4WT(African_Cubed_Chunk,[Jmax Jmax],[1 1],'forward',1);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wavelet_locs = find(Wavelet_Builder ~= 0);
% The_Big_One = zeros(2^N,2^N);
% The_Big_One(Wavelet_locs) = Orig_Wavelet_Transform(Wavelet_locs);
% subplot(3,3,6)
% GiveItBackYouMonster = angularD4WT(The_Big_One,[Jmax Jmax],[1 1],'inverse',1);
% hhhhh=imagefnan([1 1],[2^N 2^N],GiveItBackYouMonster,colmap,dax,[],[],0);
% 
% subplot(3,3,7)
% resid = GiveItBackYouMonster - Structure_of_Interest(:,:,face);
% contourf(resid)
% hold on
% plot([xbound(1) xbound(2) xbound(2) xbound(1) xbound(1)],[ybound(1) ybound(1) ybound(2) ybound(2) ybound(1)],'r')
% colorbar
% caxis([-5000 1])
% 
% subplot(3,3,8)    
% Residual_Ver2 = Maybe_Res -Structure_of_Interest(:,:,face)
% contourf(Residual_Ver2)
% hold on
% plot([xbound(1) xbound(2) xbound(2) xbound(1) xbound(1)],[ybound(1) ybound(1) ybound(2) ybound(2) ybound(1)],'r')
% colorbar
% caxis([-5000 1])