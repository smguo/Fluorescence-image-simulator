function F = imsequence_play(imsequence, varargin)  
%varagin is the timestep between each frame
max_pix_value = max(max(max(imsequence)));
min_pix_value = min(min(min(imsequence)));
map = 'gray';
ts = 0 ;
if nargin ==3, 
    ts = varargin{1} ;
    if strcmp(varargin{2},'c')
    map = 'jet';
    end
end

    if nargin == 2
       if strcmp(varargin{1},'c')
       map = 'jet';
       else          
        ts = varargin{1} ;
       end
    end


num_frame = size(imsequence,3);


for j = 1:num_frame ;
    A= imsequence(:,:,j);
    imshow(A,[],'InitialMagnification', 'fit');
% imshow(A,[]);
   
%     axis on
    colormap(map)
%     h =  colorbar ;
%     caxis([min_pix_value max_pix_value])
    
        %     imshow(A, []);
     F(j) = getframe;
title(['frame = ' num2str(j)]);
    pause(ts) ;
    drawnow
end

%  movie2avi(F,'.\movie14','FPS',10)
end