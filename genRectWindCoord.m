function rect_window = genRectWindCoord(sourceSize, sideRes, binSize)
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% This function generates a matrix rect_window which contains the
% coordinates of all the rectangular windows across an xray matrix of size
% sourceSize x sourceSize pixels that
% would be needed to scan it from the left->right and top->bottom.
% This matrix can be entered as input to the function coordClustersHist to
% generate the histogram for all windows.
% The generated windows have the form:
% The matrix is of this form:
%
%           1_r1  |  1_c1  |  1_r2  |  1_c2  |
%           2_r1  |  2_c1  |  2_r2  |  2_c2  |
%             .   |    .   |    .   |    .   |
%             .   |    .   |    .   |    .   |
%           N_r1  |  N_c1  |  N_r2  |  N_c2  |

% Where i_r1 is row and i_c1 is the column of the upper-left
% corner of a rectangle "i" that's within the dimensions of the
% xray object matrix, and i_r2 is row and i_c2 is the column of 
% the bottom-right corner of that same rectangle "i" (note that this is 
% equivalent to the field "Position" in a normal rectangle).  
% The rectangular window is generated in such a way that it scans the xray
% object matrix (which is not given as input, however, this function only
% creates the coordinates for use in function coordClustersHist) from left
% to right and top to bottom in equal horizontal and vertical steps defined
% by the variable binShift = (sourceSize-binSize)/(sideRes-1). sideRes is the number of
% bins per side that we want in the final image, and binSize is the size in
% pixels of each bin. sourceSize should be the number of rows (or columns)
% of a square matrix that we want to tile (normally this is 2048). Note
% that this function only tiles square matrices (maybe new version later).
% If binShift == binSize, then we are in the no-overlap case, and if also
% binSize*sideRes == sourceSize we are in the tesselation case (complete
% tiling with no wasted pixels and no overlap).
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
arguments
    sourceSize
    sideRes
    binSize
end
    % Define the shift between consecutive bins (both horizontally and
    % vertically)
    binShift = floor((sourceSize-binSize)/(sideRes-1));
    if (binSize*sideRes < sourceSize) 
        error("Cannot achieve total tiling of the source. \nbinSize*sideRes is smaller than the source size. \n ")
    end
    if floor((sourceSize-binSize)/(sideRes-1))~=ceil((sourceSize-binSize)/(sideRes-1)) % it means that binShift is not an integer
        fprintf("\nCannot achieve total tiling of the source without either wasting pixels or getting homogeneous resolution.\n ")
        fprintf("There are %d wasted pixels.\n",sourceSize-(binShift*(sideRes-1)+binSize)) % In this case we choose to waste some pixels on the right and bottom of the source
    end
    if mod(sourceSize,binSize)~=0 % Not super sure this is actually true
        fprintf("WARNING: The bin size is not a divisor of the source size, therefore the tiling won't be complete.\n")
    end
    % Initialize the first window 
    r1 = 1;
    c1 = 1;
    r2 = r1 + binSize -1;
    c2 = c1 + binSize -1;
    % Initialize the rect window matrix (for performance, since we altready know how big it will be anyway)
    rect_wind_temp = zeros (sideRes^2,4);
    
    % Add the first window to the matrix
    rect_wind_temp(1,:) = [r1 c1 r2 c2];
    
    % Now we start constructing the windows with a two nested for loops
    % that increment the coordinates of the rect_window
    r1_temp = r1;
    r2_temp = r2;
    for i = 1:sideRes
        for j = 1: sideRes 
            if j == 1 && i == 1
                % Add the first window to the matrix
                rect_wind_temp(1,:) = [r1 c1 r2 c2];
            else
                c1_temp = c1 + (j-1)*binShift;
                c2_temp = c2 + (j-1)*binShift;
                rect_wind_temp((i-1)*sideRes+j,:) = [r1_temp c1_temp r2_temp c2_temp]; 
            end
            
        end
        r1_temp = r1 + (i)*binShift;
        r2_temp = r2 + (i)*binShift;
    end
    rect_window = rect_wind_temp;
end