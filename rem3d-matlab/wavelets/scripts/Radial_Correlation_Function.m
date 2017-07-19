% Extends Simons' work to take apart models at different scales
% and make Radial Correlation Functions at different scales. 
% Just for fun!

% Anant Hariharan





%%%%%%%%%PREP WORK%%%%%%%%%%%%%%%%%
DepthFile = 'GN_Depths.txt';

Depths = load(DepthFile);
Depths = round(Depths/1000);

folder_interest = '/home/anant/Software/rem3d/rem3d-matlab/wavelets/codes/MIT_Wave_Results';
Corr_Matrix_Scale4 = zeros(length(Depths),length(Depths));
Corr_Matrix_Scale43 = zeros(length(Depths),length(Depths));
Corr_Matrix_Scale432 = zeros(length(Depths),length(Depths));

N = 7;
Basis = 'D4';
Jmax = 4;

All_Reconstruct_Scale4 = zeros(2^N,2^N,6,length(Depths));
All_Reconstruct_Scale43 = zeros(2^N,2^N,6,length(Depths));
All_Reconstruct_Scale432 = zeros(2^N,2^N,6,length(Depths));

DepthsStrings =[];
for i = 1:length(Depths)
    CurrString = num2str(Depths(i));
    if length(CurrString) == 1
        DepthsStrings{i} = ['000' CurrString];
    elseif length(CurrString) == 2   
        DepthsStrings{i} = ['00' CurrString];
    elseif length(CurrString) == 3   
        DepthsStrings{i} = ['0' CurrString];
    elseif length(CurrString) == 4
        DepthsStrings{i} = [CurrString];
    end
end

%%%%%%%%%%%%%%%% Begin repeated inversions and 
if exist('All_Scales_Reconstruct4.mat', 'file') == 2
for i = 1:length(Depths)
    fname1 = [folder_interest '/' 'loris5_GN_' DepthsStrings{i} '_' num2str(N) '_' num2str(Jmax) '_' Basis '_1_1.mat'];
    results1 = load(fname1);
    vw1 = results1.vw;
    
    vw1_scale4indicesonly = getkeepindex(vw1,N,Jmax,4);
    vw1_scale43indicesonly = getkeepindex(vw1,N,Jmax,[4 3]);
    vw1_scale432indicesonly = getkeepindex(vw1,N,Jmax,[4 3 2]);
    
    vw1_scale4waveletsonly = zero_wavelets(vw1,vw1_scale4indicesonly);
    vw1_scale43waveletsonly = zero_wavelets(vw1,vw1_scale43indicesonly);
    vw1_scale432waveletsonly = zero_wavelets(vw1,vw1_scale432indicesonly);
    
    vw1_recscale4 = angularD4WT(vw1_scale4waveletsonly,[Jmax Jmax],[1 1],'inverse',1);
    vw1_recscale43 = angularD4WT(vw1_scale43waveletsonly,[Jmax Jmax],[1 1],'inverse',1);
    vw1_recscale432 = angularD4WT(vw1_scale432waveletsonly,[Jmax Jmax],[1 1],'inverse',1);
        
    All_Reconstruct_Scale4(:,:,:,i) =vw1_recscale4;
    All_Reconstruct_Scale43(:,:,:,i) =vw1_recscale43;
    All_Reconstruct_Scale432(:,:,:,i) = vw1_recscale432;
    
        %     disp(i)
        %         for j = i:length(Depths)
        %             fname2 = [folder_interest '/' 'loris5_GN_' DepthsStrings{j} '_' num2str(N) '_' num2str(Jmax) '_' Basis '_1_1.mat'];
        %             results2 = load(fname2);
        %             vw2 = results2.vw;
        %             vw2_scale4indicesonly = getkeepindex(vw2,N,Jmax,4);
        %             vw2_scale43indicesonly = getkeepindex(vw2,N,Jmax,[4 3]);
        %             vw2_scale432indicesonly = getkeepindex(vw2,N,Jmax,[4 3 2]);
        %             
        %             vw2_scale4waveletsonly = zero_wavelets(vw2,vw2_scale4indicesonly);
        %             vw2_scale43waveletsonly = zero_wavelets(vw2,vw2_scale43indicesonly);
        %             vw2_scale432waveletsonly = zero_wavelets(vw2,vw2_scale432indicesonly);
        %     
        %             vw2_recscale4 = angularD4WT(vw2_scale4waveletsonly,[Jmax Jmax],[1 1],'inverse',1);
        %             vw2_recscale43 = angularD4WT(vw2_scale43waveletsonly,[Jmax Jmax],[1 1],'inverse',1);
        %             vw2_recscale432 = angularD4WT(vw2_scale432waveletsonly,[Jmax Jmax],[1 1],'inverse',1);
        %             
        %             Temp4 = corrcoef(vw1_recscale4,vw2_recscale4);
        %             Corr_Matrix_Scale4(i,j) = Temp4(1,2);
        %             Corr_Matrix_Scale4(j,i) = Corr_Matrix_Scale4(i,j);
        %             
        %             Temp43 = corrcoef(vw1_recscale43,vw2_recscale43);
        %             Corr_Matrix_Scale43(i,j) = Temp43(1,2);
        %             Corr_Matrix_Scale43(j,i) = Corr_Matrix_Scale43(i,j);
        %             
        %             Temp432 = corrcoef(vw1_recscale432,vw2_recscale432);
        %             Corr_Matrix_Scale432(i,j) = Temp432(1,2);
        %             Corr_Matrix_Scale432(j,i) = Corr_Matrix_Scale432(i,j);
        %             
        %             j
        %         end
i
end


save('All_Scales_Reconstruct4.mat','All_Reconstruct_Scale4')
save('All_Scales_Reconstruct43.mat','All_Reconstruct_Scale43')
save('All_Scales_Reconstruct432.mat','All_Reconstruct_Scale432')

else
    load('All_Scales_Reconstruct4.mat')
    load('All_Scales_Reconstruct43.mat')
    load('All_Scales_Reconstruct432.mat')
end    
%%% Actually calculate the correlation coefficients and put them in the
%%% radial correlation functions. 

for i = 1:length(Depths)
    rec1_scale4 = All_Reconstruct_Scale4(:,:,:,i);
    rec1_scale43 = All_Reconstruct_Scale43(:,:,:,i);
    rec1_scale432 = All_Reconstruct_Scale432(:,:,:,i);
    for j = i:length(Depths)
        rec2_scale4 = All_Reconstruct_Scale4(:,:,:,j);
        rec2_scale43 = All_Reconstruct_Scale43(:,:,:,j);
        rec2_scale432 = All_Reconstruct_Scale432(:,:,:,j);
        
        
        Temp4 = corrcoef(rec1_scale4,rec2_scale4);
        Corr_Matrix_Scale4(i,j) = Temp4(1,2);
        Corr_Matrix_Scale4(j,i) = Corr_Matrix_Scale4(i,j);

        Temp43 = corrcoef(rec1_scale43,rec2_scale43);
        Corr_Matrix_Scale43(i,j) = Temp43(1,2);
        Corr_Matrix_Scale43(j,i) = Corr_Matrix_Scale43(i,j);

        Temp432 = corrcoef(rec1_scale432,rec2_scale432);
        Corr_Matrix_Scale432(i,j) = Temp432(1,2);
        Corr_Matrix_Scale432(j,i) = Corr_Matrix_Scale432(i,j);
    end
    i
end







figure(1)
subplot(1,3,1)

contourf(Corr_Matrix_Scale4)
title('Correlation: Scale 4 Reconstruction')
subplot(1,3,2)
contourf(Corr_Matrix_Scale43)
title('Correlation: Scale 4,3 Reconstruction')
subplot(1,3,3)
contourf(Corr_Matrix_Scale432)
title('Correlation: Scale 4,3,2 Reconstruction')