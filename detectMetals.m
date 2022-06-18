function metal_guess = detectMetals(sum_hist_measured, hist_background, x_hist)
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% Evaluates all the histograms in the matrix sum_hist (where every row
% corresponds to the histogram of a different window, see coordClustersHist
% function to generate these) and returns a matrix in the following form:
%
%             Titanium | Copper |  Zinc  |
%               g1_Ti  |  g1_Cu |  g1_Zn |    
%               g2_Ti  |  g2_Cu |  g2_Zn | 
%                 .    |    .   |    .   | 
%                 .    |    .   |    .   | 
%               gn_Ti  |  gn_Cu |  gn_Zn | 
% Where gi_Ti, gi_Cu and gi_Zn represent the "goodness" of the guess of the
% evaluation with regards to the ith window containing titanium, copper or
% zinc respectively. A goodness of "0" means that the function has not
% found any detectable trace of for that particular metal in that
% particular window. The function works by dividing each row in
% sum_hist_measured by the array "hist_background" (which should contain
% the histogram without the ross filter), with normalization, and then by
% taking the derivative of the ratio. If the derivative goes below a
% threshold (defined below) in the 3keV-10keV range, the value of the 
% minimum is checked to see if it is within a range around the k-edge
% energy of the metals (also defined below). If this is the case, then a
% metal is said to be detected and the ratio between the minimum and the
% threshold is stored as a value of the goodness. Once all the
% histograms have been evaluated, the metal_guess matrix is returned.
% x_hist is needed because we need to know what the range of the histograms
% is in energy to be able to identify where the energies are.
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
arguments
    sum_hist_measured
    hist_background
    x_hist
end

threshold_guess = -0.009; % If the derivative goes below this value then we have detected a metal
range_metal = 300; % [eV]  The range around which we can find the minimum to say it belongs to a K-edge 
min_range = 3e3; % [eV] The minimum range we want to detect metals in
max_range = 10e3; % [eV] The maximum range we want to detect metals in
min_index = floor(((min_range-x_hist(1))*size(x_hist,2))/(x_hist(end)-x_hist(1))); % The index of the element in x_hist that corresponds to min_range
max_index = floor(((max_range-x_hist(1))*size(x_hist,2))/(x_hist(end)-x_hist(1))); % The index of the element in x_hist that corresponds to max_range

% Now, we define the ranges for detection of titanium, copper and zinc (add
% more metals if you want!)
min_Ti = 4966.4 - range_metal;
max_Ti = 4966.4 + range_metal;
min_Cu = 8978.9 - range_metal;
max_Cu = 8978.9 + range_metal;
min_Zn = 9658.6 - range_metal;
max_Zn = 9658.6 + range_metal;

% Let's find the indices of the above
min_Ti_index = floor(((min_Ti-x_hist(1))*size(x_hist,2))/(x_hist(end)-x_hist(1)));
max_Ti_index = floor(((max_Ti-x_hist(1))*size(x_hist,2))/(x_hist(end)-x_hist(1)));
min_Cu_index = floor(((min_Cu-x_hist(1))*size(x_hist,2))/(x_hist(end)-x_hist(1)));
max_Cu_index = floor(((max_Cu-x_hist(1))*size(x_hist,2))/(x_hist(end)-x_hist(1)));
min_Zn_index = floor(((min_Zn-x_hist(1))*size(x_hist,2))/(x_hist(end)-x_hist(1)));
max_Zn_index = floor(((max_Zn-x_hist(1))*size(x_hist,2))/(x_hist(end)-x_hist(1)));

% Let's smooth the background
hist_background_smooth = smoothdata(hist_background,'gaussian',60); 

% INITIALIZE MATRICES FOR PERFORMANCE
sum_hist_measured_smooth = zeros(size(sum_hist_measured,1),size(sum_hist_measured(1,:),2));
metal_tx = zeros(size(sum_hist_measured,1),size(sum_hist_measured(1,:),2));
hist_diff = zeros(size(sum_hist_measured,1),size(sum_hist_measured(1,:),2)); % The -1 is to make sure we assign to the same dimension when we differentiate... 
 
for i = 1:size(sum_hist_measured,1) 
    sum_hist_measured_smooth(i,:) = smoothdata(sum_hist_measured(i,:),'gaussian',60); % Smooth all the histograms
    metal_tx(i,:) = (sum_hist_measured_smooth(i,:)/max(sum_hist_measured_smooth(i,:)))./(hist_background_smooth/max(hist_background_smooth)); % This holds all the transmissions of the various windows
    
    diff_metal = diff(metal_tx(i,:));
    hist_diff(i,:) = [diff_metal,diff_metal(end)]; % Differentiate the measured transmission to see where the K-edge is
    %hist_diff(i,end+1)=hist_diff(i,end); % Differentiation loses an element, so we duplicate the final one
    
    % Now that we have all the derivatives, we want to see (within the
    % range [min_range max_range]) what elements of each derivative go
    % "below" threshold (they need to be more negative then this threshold).
    below_threshold(i,:) = zeros(size(x_hist));
    % To do this, we first extract the elements below threshold from each
    % derivative in the range defined by [min_range max_range]:
    % NOTE: THE VECTORS BELOW HAVE A SMALLER DIMENSION THAN HIST_DIFF
    % Zero pad so it's the same dimension as hist_diff
    below_threshold(i,1:(min_index-1)) = zeros(1,min_index-1);
    below_threshold(i,min_index:max_index) = hist_diff(i,min_index:max_index).*(hist_diff(i,min_index:max_index)<threshold_guess);
    below_threshold(i,max_index:end) = zeros(1,size(x_hist,2)-max_index+1);
    
    %below_threshold = zeros(size(x_hist)) + below_threshold_temp; 
%     if sum(below_threshold(i,:),'all') == 0 % There were no elements of hist_diff(i) in the range which were below threshold
%         metal_guess_temp(i,:) = [0 0 0]; % There are no metals in this region, this is just background
%         break
%     end
    
    % Now we check the metal windows!
    metal_guess_temp(i,:) = [0 0 0]; % We must intialize this to increment it later
    % TITANIUM WINDOW RANGE:
    if sum(below_threshold(i,min_Ti_index:max_Ti_index))~=0
        g_Ti = min(below_threshold(i,min_Ti_index:max_Ti_index))/threshold_guess;
        metal_guess_temp(i,:) = metal_guess_temp(i,:) + [g_Ti-1 0 0];   
    end
    % COPPER WINDOW RANGE:
    if sum(below_threshold(i,min_Cu_index:max_Cu_index))~=0
        g_Cu = min(below_threshold(i,min_Cu_index:max_Cu_index))/threshold_guess;
        metal_guess_temp(i,:) = metal_guess_temp(i,:) + [0 g_Cu-1 0];   
    end
    % ZINC WINDOW RANGE:
    if sum(below_threshold(i,min_Zn_index:max_Zn_index))~=0
        g_Zn = min(below_threshold(i,min_Zn_index:max_Zn_index))/threshold_guess;
        metal_guess_temp(i,:) = metal_guess_temp(i,:) + [0 0 g_Zn-1];   
    end    
   
end
% Now we have all the guesses, but we still need to normalize them!
%metal_guess_norm
metal_guess = metal_guess_temp;



end