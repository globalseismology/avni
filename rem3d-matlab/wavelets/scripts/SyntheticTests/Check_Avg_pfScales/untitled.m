%%%%% Compare Averages of Wavelets

N = 7;
Jmax = [4 5 6 7];

[vwlev,vwlevs]=cube2scale(N,[Jmax Jmax],1);

%%%%%%

for i = 1:length(JMax)
    TransformMe = zeros(2^N,2^N);
    scalemax_indices = find(vwlev == Jmax(i));
    TransforMe(scalemax_indices) = 1;
    reconstructed_synthDv = angularD4WT(TransformMe,[Jmax Jmax],[1 1],'inverse',1);
    plotoncube(reconstructed_synthDV,'2D')
end


