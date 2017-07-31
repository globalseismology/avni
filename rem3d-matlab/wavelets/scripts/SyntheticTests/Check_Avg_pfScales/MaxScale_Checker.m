%%%%% Compare Averages of Wavelets
%Enter the list of scales you're interested in. 
N = [9 9 9 9 9 9];
Jmax = [6 6 6 6 6 6];
Scale_Interest = [1 2 3 4 5 6];
Mean_vec = [];

%%%%%%

for i = 1:length(Jmax)
    [vwlev,vwlevs]=cube2scale(N(i),[Jmax(i) Jmax(i)],1);
    TransformMe = zeros(2^N(i),2^N(i));
    scale_indices = find(vwlev == Scale_Interest(i));
    TransformMe(scale_indices) = 10;
    reconstructed_synthDv = angularD4WT(TransformMe,[Jmax(i) Jmax(i)],[1 1],'inverse',1);
    plotoncube(reconstructed_synthDv,'2D')
    CurrMean = mean(mean(reconstructed_synthDv));
    Mean_vec =[Mean_vec CurrMean];
end


