function metalChannels = genMultispectralImage(metal_guess,rect_windows)
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% This function will generate as many images as the columns in the matrix
% metal_guess, which is the output of detectMetals function. In a first
% example, detectMetals can detect Titanium, Copper and Zinc K-edges, so
% metal_guess has 3 columns and metal_channels has 3 elements. The ith 
% element of metalChannels is the image corresponding to the ith metal in
% metal_guess, ordered by increasing K-edge. Therefore, metalChannels is a
% 3D matrix where metalChannels(1,:,:) is the Titanium channel, 
% metalChannels(2,:,:) is the Copper channel and metalChannels(3,:,:) is 
% the Zinc channel.
% The size of each channel image is defined by the square root of the size
% of rect_windows (a property called sideRes in the input of
% genRectWindCoord defines this quantity). 
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
arguments
    metal_guess
    rect_windows
end
sideRes = sqrt(size(rect_windows,1));
metalChannels_temp = zeros(sideRes,sideRes,3); % Initialize the metalChannels matrix

for row = 1:sideRes 
    for col = 1:sideRes 
        % TITANIUM CHANNEL
        metalChannels_temp(row,col,1) = metal_guess((row-1)*sideRes+col,1);
        % COPPER CHANNEL
        metalChannels_temp(row,col,2) = metal_guess((row-1)*sideRes+col,2);
        % ZINC CHANNEL
        metalChannels_temp(row,col,3) = metal_guess((row-1)*sideRes+col,3);
        
    end
end
metalChannels = metalChannels_temp;
end