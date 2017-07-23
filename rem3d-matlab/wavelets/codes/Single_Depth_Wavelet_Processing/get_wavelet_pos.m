function [ wavelet_locs ] = get_wavelet_pos( x,y,N,Jmax,face,plot  )
% Get array of points with binary signatures 
% where a wavelet exists, given a point on the cubed sphere. 

% Note that we use rounding to ensure that we actually generate indices. 
defval('plot',0)

orig_locs = zeros(2^N,2^N,6);
for i = 1:length(x)
orig_locs(x(i),y(i),face)= 1;
end
wavelet_locs= zeros(2^N,2^N,6);

Length = 2^N;

for i = 1:Jmax
    for j = 1:length(x)
        xpts_corner1 = round(x(j)/(2^i) + Length/(2^i));
        ypts_corner1 =  round(y(j)/(2^i));

        xpts_corner1(xpts_corner1 == 0) = 1;
        ypts_corner1(ypts_corner1 == 0) = 1;


        xpts_corner2 = round(x(j)/(2^i));
        ypts_corner2 = round(y(j)/(2^i) + Length/(2^i));


        xpts_corner2(xpts_corner2 == 0) = 1;
        ypts_corner2(ypts_corner2 == 0) = 1 ;  


        xpts_corner3 =round( x(j)/(2^i) + Length/(2^i));
        ypts_corner3 = round(y(j)/(2^i) + Length/(2^i));


        xpts_corner3(xpts_corner3 == 0) = 1;
        ypts_corner3(ypts_corner3 == 0) = 1; 

        wavelet_locs(xpts_corner1,ypts_corner1,face) = 1;
        %%%
        if xpts_corner1-1 ~= 0 && ypts_corner1-1 ~= 0 && xpts_corner1-1 <= 2^N && ypts_corner1-1 <= 2^N
            wavelet_locs(xpts_corner1-1,ypts_corner1-1,face) = 1;    
        end
        
        if xpts_corner1-1 ~= 0 && ypts_corner1 ~= 0 && xpts_corner1-1 <= 2^N && ypts_corner1 <= 2^N
            wavelet_locs(xpts_corner1-1,ypts_corner1,face) = 1;    
        end
        
        if xpts_corner1-1 ~= 0 && ypts_corner1+1 ~= 0 && xpts_corner1-1 <= 2^N && ypts_corner1+1 <= 2^N
            wavelet_locs(xpts_corner1-1,ypts_corner1+1,face) = 1;    
        end
        
        %%%
        if xpts_corner1 ~= 0 && ypts_corner1-1 ~= 0 && xpts_corner1 <= 2^N && ypts_corner1-1 <= 2^N
            wavelet_locs(xpts_corner1,ypts_corner1-1,face) = 1;    
        end
        
        if xpts_corner1 ~= 0 && ypts_corner1+1 ~= 0 && xpts_corner1 <= 2^N && ypts_corner1+1 <= 2^N
            wavelet_locs(xpts_corner1,ypts_corner1+1,face) = 1;    
        end
        
        %%%
        if xpts_corner1+1 ~= 0 && ypts_corner1-1 ~= 0 && xpts_corner1+1 <= 2^N && ypts_corner1-1 <= 2^N
            wavelet_locs(xpts_corner1+1,ypts_corner1-1,face) = 1;    
        end
        
        if xpts_corner1+1 ~= 0 && ypts_corner1 ~= 0 && xpts_corner1+1 <= 2^N && ypts_corner1 <= 2^N
            wavelet_locs(xpts_corner1+1,ypts_corner1,face) = 1;    
        end
            
        if xpts_corner1+1 ~= 0 && ypts_corner1+1 ~= 0 && xpts_corner1+1 <= 2^N && ypts_corner1+1 <= 2^N
            wavelet_locs(xpts_corner1+1,ypts_corner1+1,face) = 1;    
        end
        %%%
        
        
        
        wavelet_locs(xpts_corner2,ypts_corner2,face) = 1;
        
        
        %%%
        if xpts_corner2-1 ~= 0 && ypts_corner2-1 ~= 0 && xpts_corner2-1 <= 2^N && ypts_corner2-1 <= 2^N
            wavelet_locs(xpts_corner2-1,ypts_corner2-1,face) = 1;    
        end
        
        if xpts_corner2-1 ~= 0 && ypts_corner2 ~= 0 && xpts_corner2-1 <= 2^N && ypts_corner2 <= 2^N
            wavelet_locs(xpts_corner2-1,ypts_corner2,face) = 1;    
        end
        
        if xpts_corner2-1 ~= 0 && ypts_corner2+1 ~= 0 && xpts_corner2-1 <= 2^N && ypts_corner2+1 <= 2^N
            wavelet_locs(xpts_corner2-1,ypts_corner2+1,face) = 1;    
        end
        
        %%%
        if xpts_corner2 ~= 0 && ypts_corner2-1 ~= 0 && xpts_corner2 <= 2^N && ypts_corner2-1 <= 2^N
            wavelet_locs(xpts_corner2,ypts_corner2-1,face) = 1;    
        end
        
        if xpts_corner2 ~= 0 && ypts_corner2+1 ~= 0 && xpts_corner2 <= 2^N && ypts_corner2+1 <= 2^N
            wavelet_locs(xpts_corner2,ypts_corner2+1,face) = 1;    
        end
        
        %%%
        if xpts_corner2+1 ~= 0 && ypts_corner2-1 ~= 0 && xpts_corner2+1 <= 2^N && ypts_corner2-1 <= 2^N
            wavelet_locs(xpts_corner2+1,ypts_corner2-1,face) = 1;    
        end
        
        if xpts_corner2+1 ~= 0 && ypts_corner2 ~= 0 && xpts_corner2+1 <= 2^N && ypts_corner2 <= 2^N
            wavelet_locs(xpts_corner2+1,ypts_corner2,face) = 1;    
        end
            
        if xpts_corner2+1 ~= 0 && ypts_corner2+1 ~= 0 && xpts_corner2+1 <= 2^N && ypts_corner2+1 <= 2^N
            wavelet_locs(xpts_corner2+1,ypts_corner2+1,face) = 1;    
        end
        %%%
        
        
        
        wavelet_locs(xpts_corner3,ypts_corner3,face) = 1;
        
        
        %%%
        if xpts_corner3-1 ~= 0 && ypts_corner3-1 ~= 0 && xpts_corner3-1 <= 2^N && ypts_corner3-1 <= 2^N
            wavelet_locs(xpts_corner3-1,ypts_corner3-1,face) = 1;    
        end
        
        if xpts_corner3-1 ~= 0 && ypts_corner3 ~= 0 && xpts_corner3-1 <= 2^N && ypts_corner3 <= 2^N
            wavelet_locs(xpts_corner3-1,ypts_corner3,face) = 1;    
        end
        
        if xpts_corner3-1 ~= 0 && ypts_corner3+1 ~= 0 && xpts_corner3-1 <= 2^N && ypts_corner3+1 <= 2^N
            wavelet_locs(xpts_corner3-1,ypts_corner3+1,face) = 1;    
        end
        
        %%%
        if xpts_corner3 ~= 0 && ypts_corner3-1 ~= 0 && xpts_corner3 <= 2^N && ypts_corner3-1 <= 2^N
            wavelet_locs(xpts_corner3,ypts_corner3-1,face) = 1;    
        end
        
        if xpts_corner3 ~= 0 && ypts_corner3+1 ~= 0 && xpts_corner3 <= 2^N && ypts_corner3+1 <= 2^N
            wavelet_locs(xpts_corner3,ypts_corner3+1,face) = 1;    
        end
        
        %%%
        if xpts_corner3+1 ~= 0 && ypts_corner3-1 ~= 0 && xpts_corner3+1 <= 2^N && ypts_corner3-1 <= 2^N
            wavelet_locs(xpts_corner3+1,ypts_corner3-1,face) = 1;    
        end
        
        if xpts_corner3+1 ~= 0 && ypts_corner3 ~= 0 && xpts_corner3+1 <= 2^N && ypts_corner3 <= 2^N
            wavelet_locs(xpts_corner3+1,ypts_corner3,face) = 1;    
        end
            
        if xpts_corner3+1 ~= 0 && ypts_corner3+1 ~= 0 && xpts_corner3+1 <= 2^N && ypts_corner3+1 <= 2^N
            wavelet_locs(xpts_corner3+1,ypts_corner3+1,face) = 1;    
        end
        %%%
                
    end
end
    
    

for k = 1:length(x)
    %scaling_function
    xpts_corner = round(x(k)/(2^Jmax)); 
    ypts_corner = round(y(k)/(2^Jmax));
    xpts_corner(xpts_corner == 0) = 1;
    ypts_corner(ypts_corner == 0) = 1; 

    wavelet_locs(xpts_corner,ypts_corner,face) = 1;
    
    
    %%%
    if xpts_corner+1 ~= 0 && ypts_corner-1 ~= 0 && xpts_corner+1 <= 2^N && ypts_corner-1 <= 2^N
        wavelet_locs(xpts_corner+1,ypts_corner-1,face) = 1;    
    end
        
    if xpts_corner+1 ~= 0 && ypts_corner ~= 0 && xpts_corner+1 <= 2^N && ypts_corner <= 2^N
        wavelet_locs(xpts_corner+1,ypts_corner,face) = 1;    
    end
            
    if xpts_corner+1 ~= 0 && ypts_corner+1 ~= 0 && xpts_corner+1 <= 2^N && ypts_corner+1 <= 2^N
        wavelet_locs(xpts_corner+1,ypts_corner+1,face) = 1;    
    end
        %%%
    
    
    
    
    
end

if plot ==1
subplot(1,2,1);    
h=imagefnan([1 1],[2^N 2^N],wavelet_locs(:,:,face));
subplot(1,2,2);
g=imagefnan([1 1],[2^N 2^N],orig_locs(:,:,face));
end

end