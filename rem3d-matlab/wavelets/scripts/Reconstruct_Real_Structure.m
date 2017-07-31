
clear all
clc
Look_At_Residual = 0;
%Non-Synthetic Tests: Try on Real structure. For this script, african
%topography! 
compare_structure = 1;
setenv('IFILES','/home/moulik/Software/fjsimons-MMASS-v1.0.1/DATA')
colperc = [5 95];
Get_Edgy = 0;
N = 7; Jmax = 4; face = 3; J = 4;
defval('L',ceil(2^(N+1)))
defval('colmap','kelicol')
xbound =[40 90];
ybound = [60 120];
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

% %heh
% v(:,:,face) = 1;
% African_Cubed_Chunk(:,:) = 1;
% 
% %




dax=prctile(African_Cubed_Chunk(:),colperc);
subplot(3,3,1)
h=imagefnan([1 1],[2^N 2^N],v(:,:,face),colmap,dax,[],[],0);
title('original anomaly')

vw=angularD4WT(v,[Jmax Jmax],[1 1],'forward',1);
subplot(3,3,2)
i=imagefnan([1 1],[2^N 2^N],vw(:,:,face))

title('wavelet coefficients describing anomaly')

Reconstruct_Africa_Direct = angularD4WT(vw,[J J],[1 1],'inverse',1);
subplot(3,3,3)
j=imagefnan([1 1],[2^N 2^N],Reconstruct_Africa_Direct(:,:,face),colmap,dax)
%xlim(xbound)
%ylim(ybound)
title('Immediate Reconstruction')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Let's say the EARS is around x = 40:90, y = 60:120 on cubed sphere for N =
%7. Hard thresholding. 
[ xpts,ypts ] = Get_Interior_PolygonPts_on_CubedSphere( N,[ybound(1) ybound(1) ybound(2) ybound(2)],[xbound(1) xbound(2) xbound(2) xbound(1)] );
subplot(3,3,[4,5])
plot(xpts,ypts,'ro')
title('points of interest')


%Detect Wavelets where we want them to be. 
if Get_Edgy == 1;
    [ wavelet_locs ] = get_wavelet_pos( xpts,ypts,N,Jmax,face,0,1);
else
    [ wavelet_locs ] = get_wavelet_pos( xpts,ypts,N,Jmax,face,0,0);

end

subplot(3,3,6)
k=imagefnan([1 1],[2^N 2^N],wavelet_locs(:,:,face));
title('Location of Wavelet Coefficients correponding to region of interest')

%Now, get those wavelets~

[wavelets_of_interest ] = zero_wavelets(vw,wavelet_locs,'db4',N,1,Jmax );
[maybe_Result] = angularD4WT(wavelets_of_interest,[Jmax Jmax],[1 1],'inverse',1);

subplot(3,3,7)

l=imagefnan([1 1],[2^N 2^N],wavelets_of_interest(:,:,face));
title('Zeroed Out Wavelets')

subplot(3,3,7)

l=imagefnan([1 1],[2^N 2^N],wavelets_of_interest(:,:,face));
title('Zeroed Out Wavelets')

subplot(3,3,9)

m=imagefnan([1 1],[2^N 2^N],maybe_Result(:,:,face),colmap,dax,[],[],0);
%xlim(xbound)
%ylim(ybound)
title('Reconstructed Geographical Region of Interest')








%%%%%%%%%%%%%%%%%% Compare Structure
if compare_structure ==1
    
    Structure_of_Interest = zeros(2^N,2^N,6);
    for i = 1:length(xpts)
        
        Structure_of_Interest(xpts(i),ypts(i),face) = Structure_of_Interest(xpts(i),ypts(i),face) + Reconstruct_Africa_Direct(xpts(i),ypts(i),face);
    end
    figure
    subplot(2,2,1)
    aa=imagefnan([1 1],[2^N 2^N],Structure_of_Interest(:,:,face),colmap,dax,[],[],0);

    title('Truncated, Original Geographical Region of Interest')
    
    subplot(2,2,2)

    m=imagefnan([1 1],[2^N 2^N],maybe_Result(:,:,face),colmap,dax,[],[],0);
    %xlim(xbound)
    %ylim(ybound)
    title('Reconstructed Geographical Region of Interest')




    
    %%Make Difference. 

    Difference = maybe_Result(:,:,face) - Structure_of_Interest(:,:,face);
    subplot(2,2,3)
    %bb=imagefnan([1 1],[2^N 2^N],Difference,colmap,dax,[],[],0);
    contourf(Difference) 
    colorbar
    hold on
    plot([xbound(1) xbound(2) xbound(2) xbound(1) xbound(1)],[ybound(1) ybound(1) ybound(2) ybound(2) ybound(1)],'r')
    title('Residual of Geographic Reconstruction wrt Original')
    
    Norm = Difference.*Difference;
    subplot(2,2,4)
    %cc=imagefnan([1 1],[2^N 2^N],Norm);
    contourf(Norm)
    colorbar
    hold on
    plot([xbound(1) xbound(2) xbound(2) xbound(1) xbound(1)],[ybound(1) ybound(1) ybound(2) ybound(2) ybound(1)],'r')
    title('L2 (sorta) Norm of Geographic Reconstruction wrt Original')

    
    

end
%%%%%%%%%%

%Analyse difference
Residual_WaveletDecomp = angularD4WT(Difference,[Jmax Jmax],[1 1],'forward',1);
figure
dd=imagefnan([1 1],[2^N 2^N],Residual_WaveletDecomp,colmap,dax,[],[],0);
colorbar
caxis([-1 1])
hold on


if Look_At_Residual ==1
    figure
for i = 1:2^N
    for j = 1:2^N
        if wavelet_locs(i,j,face) ~= 0
            scatter(i,j,'r+')
            hold on
        end
    end
end

for i = 1:2^N
    for j = 1:2^N
        if Residual_WaveletDecomp(i,j) ~= 0
            scatter(i,j,5,Residual_WaveletDecomp(i,j))
            hold on
        end
    end
end
end
%%%%%%%%%%%

%%compare inverted truncation....?

  figure
  subplot(1,3,1)
  aaa=imagefnan([1 1],[2^N 2^N],Structure_of_Interest(:,:,face),colmap,dax,[],[],0);
  title('Truncated, Original Geographical Region of Interest')
  
  Wavelet_Coeffs_Truncated_Original = angularD4WT(Structure_of_Interest,[J J],[1 1],'forward',1);
  subplot(1,3,2)
  bbb=imagefnan([1 1],[2^N 2^N],Wavelet_Coeffs_Truncated_Original(:,:,face),colmap,dax,[],[],0);
  title('Wavelet coeffs of Truncated region')
  
  Truncated_Reconstruction = angularD4WT(Wavelet_Coeffs_Truncated_Original,[J J],[1 1],'inverse',1);
  subplot(1,3,3)
  ccc=imagefnan([1 1],[2^N 2^N],Truncated_Reconstruction(:,:,face),colmap,dax,[],[],0);
  title('Wavelet coeffs of Truncated region')
  
%  %%%%%%%%%%%%%%%%%%%%%%%
%  figure
%  for i = 1:2^N
%     for j = 1:2^N
%         if wavelet_locs(i,j,face) ~= 0
%             scatter(i,j,'r+')
%             hold on
%         end
%     end
% end
% 
% for i = 1:2^N
%     for j = 1:2^N
%         if Wavelet_Coeffs_Truncated_Original(i,j,face) ~= 0
%             scatter(i,j,5,Wavelet_Coeffs_Truncated_Original(i,j,face),'s')
%             hold on
%         end
%     end
% end
%  

%%%%% Locate with Truncated Coeffs
figure
  Wavelet_Coeffs_Truncated_Original = angularD4WT(Structure_of_Interest,[J J],[1 1],'forward',1);
  subplot(2,3,1)
  bbb=imagefnan([1 1],[2^N 2^N],Wavelet_Coeffs_Truncated_Original(:,:,face),colmap,dax,[],[],0);
  title('Wavelet coeffs of Truncated region')
  
  Truncated_Struc_vwloc = zeros(128,128,6);
  a = find(Wavelet_Coeffs_Truncated_Original ~= 0);
  Truncated_Struc_vwloc(a) = 1;
  subplot(2,3,2)
  for i = 1:2^N
      for j = 1:2^N
          if Truncated_Struc_vwloc(i,j,face) ==1
          scatter(i,j,Truncated_Struc_vwloc(i,j,face))
          end
            hold on
      end
  end
          
[wavelets_of_interest_try2 ] = zero_wavelets(vw,Truncated_Struc_vwloc,'db4',N,1,Jmax );
 subplot(2,3,4)
 bbb=imagefnan([1 1],[2^N 2^N],wavelets_of_interest_try2(:,:,face),colmap,dax,[],[],0);
 title('Actual Value: Wavelet coeffs of Truncated region')
 
 
out_v = angularD4WT(wavelets_of_interest_try2,[J J],[1 1],'inverse',1);
  subplot(2,3,3)
  ddd=imagefnan([1 1],[2^N 2^N],out_v(:,:,face),colmap,dax,[],[],0);
  title('Reconstructed Truncated region using locations from inverted structure')
  
  
  
subplot(2,3,[5,6])  
Difference3 = out_v(:,:,face) - Structure_of_Interest(:,:,face);

hold on
contourf(Difference3)
plot([xbound(1) xbound(2) xbound(2) xbound(1) xbound(1)],[ybound(1) ybound(1) ybound(2) ybound(2) ybound(1)],'r')
    colorbar
    hold on
    
