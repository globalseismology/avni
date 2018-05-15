%%%SETTLE THIS ISSUE! HOW TO RECONCILE ORIG-VS-WAVELET STRUCTURE?

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

vw=angularD4WT(v(:,:,face),[Jmax Jmax],[1 1],'forward',1);
ReturnAfricaToMe = angularD4WT(vw,[Jmax Jmax],[1 1],'inverse',1);
dax=prctile(African_Cubed_Chunk(:),colperc);
subplot(2,3,1)
h=imagefnan([1 1],[2^N 2^N],African_Cubed_Chunk,colmap,dax,[],[],0);
title('Orig Structure')
subplot(2,3,2)
h=imagefnan([1 1],[2^N 2^N],ReturnAfricaToMe,colmap,dax,[],[],0);
title('Recovered Structure- inv back and forth')
Residual = African_Cubed_Chunk - ReturnAfricaToMe;
subplot(2,3,3)
contourf(Residual)
colorbar
colormap(redbluecmap);
caxis([-0.000000001 0.000000001])
title('Residual of Original and wv reconstruction')

Comp1 = out_v(:,:,face) -African_Cubed_Chunk;
Comp2 = out_v(:,:,face) - ReturnAfricaToMe;

subplot(2,3,4)
contourf(Comp1)
% for i = 1:2^N
%     for j = 1:2^N
%         scatter(i,j,5,100*Comp1(i,j))
%         hold on
%     end
% end
title('Residual of EARS wrt Original Structure')
colorbar
caxis([-20 20])
%caxis([-1000 1000])
subplot(2,3,5)
contourf(Comp2)
% for i = 1:2^N
%     for j = 1:2^N
%         scatter(i,j,5,100*Comp2(i,j))
%         hold on
%     end
% end
title('Residual of EARS wrt wavelet reconstruction')
colorbar
caxis([-20 20])
%caxis([-100 100])

%%%%%Find max residual
maxresid = 0;
for i = 60:120
    for j = 40:90
        resid = Comp1(i,j);
        if abs(resid)>maxresid
            maxresid = abs(resid);
        end
    end
end
disp(maxresid)