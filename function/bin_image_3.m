% no interpolation. the residue are cut off
function imBinned = bin_image_3(im, binning)
    Lx = size(im,1);
    Ly = size(im,2);
    Lz = size(im,3);
%     if mod(Lx,binning)==0 && mod(Ly,binning)==0 % true binning
        
        imBinned = zeros(floor(Lx/binning),floor(Ly/binning), Lz);
        for ix = 1 : size(imBinned,1)
            for iy = 1: size(imBinned,2)
                xIndexList = (binning*(ix-1)+1) : (ix*binning);
                yIndexList = (binning*(iy-1)+1) : (iy*binning);
                imBinned(ix,iy,:)=sum(sum(im(xIndexList,yIndexList,:),1),2);
            end
        end
%     else % linear interpolation
%         imBinned = imresize(im,1/binning,'bilinear')*binning*binning;

%     end

end

