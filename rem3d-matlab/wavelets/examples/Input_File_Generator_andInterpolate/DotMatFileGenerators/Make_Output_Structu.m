%%%%

%Make a map of dv or inverse anomalies. 
N = 7 
Jmax = 4;
eo = 1;
face = 1;
inv = 0;


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
        if inv == 1
        dv_map = sparse((angularD4WT(pts(:,:,face),[Jmax Jmax],[1 1],'forward',1)));
        elseif inv == 0
        dv_map = sparse((angularD4WT(pts(:,:,face),[Jmax Jmax],[1 1],'inverse',1)));
        end
        
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

for i = 1:2^N
    for j = 1:2^N
        linearInd = sub2ind([2^N 2^N],i,j);
        Out_Struc(linearInd).dv = All_map(i).map{j};
    end
end
Out_Struc(1).MetaN = N;
Out_Struc(1).MetaJmax = Jmax;
Out_Struc(1).Metaeo = eo;


if inv == 1;
Out_Struc(1).MetaForwardBackward = 'WCoeffs';
elseif inv ==0;
Out_Struc(1).MetaForwardBackward = 'Mmaps';
end   
Container.Me = Out_Struc;
save([Out_Struc.MetaForwardBackward '.N' num2str(N) '.J' num2str(Jmax) '.D4.mat'],'-struct','Container')



%%%%%%%Reformat the array to make it easier to work with if necessary, but
%%%%%%%so far...
