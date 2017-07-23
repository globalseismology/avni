%%%%

%Make a map of dv anomalies. 
N = 7 
Jmax = 4;
face = 1;



%Things get weird because of parfor: can only create variables within the
%for loops, forcing us to re-initialize the pts variable every loop :(



for i = 1:2^N
   
    [vwlev,vwlevs]=cube2scale(N,[Jmax Jmax],1);
    
    %Make the array which we will store results in- parfor likes the cell
    %arrays, (according to the wise people on the internet!- I should check this!) 
    mat_array = cell(1,2^N);
    parfor j = 1:2^N
        counter = (i-1)*(2^N)+j;
        pts = zeros(2^N,2^N,6);
        pts(i,j,face) = 1;
        % Make the results of inversion a sparse matrix. Invert a single
        % point, defined by i,j. 
        dv_map = sparse((angularD4WT(pts(:,:,face),[Jmax Jmax],[1 1],'inverse',1)));
        
        %Put each one in the cell array.
        mat_array{1,j} =dv_map;
        
        %Reset pts, but this really doesn't matter since we hard reset it
        %at start.
        pts(i,j,face) = 0;
        
        %Visually see progress on command line. 
        disp([i j])
    end
    All_map(i).map = mat_array;
end

save(['invwavelet.' num2str(N) '.' num2str(Jmax) '.D4.mat'],'All_map')



%%%%%%%Reformat the array to make it easier to work with if necessary, but
%%%%%%%so far...

