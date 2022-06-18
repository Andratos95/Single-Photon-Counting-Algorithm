classdef xrayMatrix < handle
    properties
        A % The image matrix (unfiltered)
        info % Informations on the image
        numRows % The image matrix height in pixels
        numCols % The image matrix width in pixels
        visMatrix % The visited elements matrix
        clusterMatrix % The cluster matrix
        noiseVal % The average amount of noise in the image
        thresholdVal % The value of the threshold to discriminate clusters from noise
        fileName % The name of the file as a string (set to 0 if just a matrix)
        coordClusters % The matrix of coordinates of the clusters (initialized to [0 0 0 0], find with findCoordClusters)
        
    end
    methods (Static)
        function y = readXray(im_filename) 
            % This static method is a utility function in case one wants to
            % transorm a .tif image into a matrix (not so useful actually
            % after the xrayMatrix v2 has been implemented which can accept
            % both types of inputs...
            arguments
                im_filename % This should be a string, the name of the image you want to read
            end
            y = double(imread(im_filename)); % Converts image to matrix of doubles
            
        end
    end
    methods
        function [obj] = xrayMatrix(A) % THIS IS AN OBJECT CONSTRUCTOR (A is an image matrix OR a .tif image filename)
            if isa(A,'string')||isa(A,'char') % You can choose to either input a filename (a string or sequence of char) 
                obj.A = double(imread(A));
                obj.info = imfinfo(A); % Retrieves metadata information from the image
                obj.numRows = obj.info.Height; % Updates numRows to the height of the image (# of rows)
                obj.numCols = obj.info.Width; % Updates numCols to the width of the image (# of columns)
                
                obj.fileName = A; % This just keeps the name of the file as part of the object so the methods can know what it is
                                     % and use it as necessary.
            elseif isa(A,"double") % But you also can just input its matrix (this is useful for submatrices)
                obj.A = A; % Converts image to matrix of doubles
                obj.numRows = size(A,1); % Updates numRows to the height of the image (# of rows)
                obj.numCols = size(A,2); % Updates numCols to the width of the image (# of columns)
                obj.info = 0; % If this field is 0, it means the xrayMatrix object has no info (it does not derive from a .tif image...)
                obj.fileName = 0; % If this field is 0, it means the xrayMatrix object has no correspondent filename
                                  
            else
                msg = "You must either insert an image (as a string) or a matrix. The other field should be set to 0.";
                error(msg)
            end
            
            obj.visMatrix = zeros(obj.numRows,obj.numCols); % Matrix of the visited elements. It has the same dimension as A, and
                                                            % is initialized with zeros. If the elementh ij becomes 1,
                                                            % it means that the ij element of A has been visited by
                                                            % the search algorithm.
            obj.clusterMatrix = zeros(obj.numRows,obj.numCols); % Matrix that keeps track if a cell is part of a cluster,
                                                                % used to avoid counting elements as clusters more than once.
            obj.noiseVal = 301; % This is the electronic offset bias in the photo
            obj.thresholdVal = 320; % Set this threshold some standard deviation above offset (read-out noise around 2.5 counts @-70°)
            obj.coordClusters = [0 0 0 0];
        end
        function [] = displayXray(obj,lowLim,highLim) 
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% Displays the elements of the image above threshold and below upper limit.
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
            arguments
                obj 
                lowLim (1,1) 
                highLim (1,1)  
            end
            fprintf("Displaying X-ray image %s with values in range [%d, %d].\n",obj.fileName,lowLim,highLim);
            figure()  
            imshow(obj.A.*(obj.A>lowLim & obj.A<highLim),[]);
            colormap(hot)
        end 
        function [] = displayColumn(obj,columnNumber)
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% Displays a 2D plot of the matrix along a column defined by columnNumber
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
            arguments
                obj 
                columnNumber (1,1) 
            end
            fprintf("Plotting column number %d\n",columnNumber);
            figure 
            plot(1:obj.info.Width,obj.A(:,columnNumber)) % Plots a column of the matrix (use to visually estimate noise level) 
            axis equal
            title("Line cut plot","FontSize",20)
            txt = ["Column: " num2str(columnNumber)];
            subtitle(txt)
            xlabel("pixels","FontSize",12)
            ylabel("counts","FontSize",12)
        end
        function sat_image = saturateXrays(obj)
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% This function will sum all the matrices of the vector obj and
% will return the saturated image as an xray object. 
% Make sure all images have the same dimensions!!!
% You can use displayXray function to display the image.
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
            arguments
                obj
            end
            M = zeros(obj(1).numRows,obj(1).numCols);
            for i = 1:size(obj,2)
                M = M + obj(i).A;
            end
            sat_image = xrayMatrix(M);
        end
        function y = isVisited(obj,i,j) 
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% Checks whether the object's cell of coordinates (i,j) has
% been previously visited. Returns a boolean ("true" if
% visited, "false" if unvisited).
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
            arguments
                obj 
                i (1,1) 
                j (1,1) 
            end
            if obj.visMatrix(i,j)==1
                y=true;
            else
                y=false;
            end
        end
        function [] = visit(obj,i,j)
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% Visits the element (i,j) of the object's visited matrix.
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
            arguments
                obj 
                i (1,1) 
                j (1,1) 
            end
            obj.visMatrix(i,j)=1;
        end
        function y = isSafe(obj, row, col)
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% Checks whether the coordinates (row,col) are within the bounds of the
% given matrix and unvisited. Returns a boolean.
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
            arguments
                obj %{mustBeA(obj,"xrayMatrix")}
                row (1,1) 
                col (1,1) 
            end

            if row <= 0 || col <= 0 || row > obj.numCols || col > obj.numRows || isVisited(obj,row,col) == true
                y = false;
            else 
                y = true;
            end
        end
        function y = isInBounds(obj, row, col)
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% Checks whether the coordinates (row,col) are within the bounds of the
% given matrix (doesn"t check for being visited). Returns a boolean 
% ("true" if in bounds, "false" if out of bounds).
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
            arguments
                obj 
                row (1,1) 
                col (1,1) 
            end
            if (row <= 0 || col <= 0 || row > obj.numCols || col > obj.numRows) == true
                y = false;
            else 
                y = true;
            end
        end
        function timeLoop (obj,row,col)
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% A utility function used to time the process of looping
% through the matrix. It gives both elapsed and an
% approximation of the remaining time. (WARNING: SLOWS DOWN
% THE CODE, therefore it should probably not be used. Also, it
% has troubles with parallel for loops).
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
            arguments
                obj 
                row (1,1)
                col (1,1)
            end
            if mod((row*10000)/(obj.numCols-mod(obj.numCols,1000)),100) == 0 && col == 1 % This occurs every approx 1% progress increment
                t = toc; % Keeps track of elapsed time
                progress_percentage = floor(100*(row/(obj.numCols-mod(obj.numCols,1000))));
                t_ave = t/progress_percentage; % Average time for 1% progress 
                fprintf("Progress: %d %% // ",progress_percentage);
                fprintf("Elapsed: %d min, %f s // Remaining: %d min, %f s.\n", floor(t/60),mod(t,60),floor((t_ave*(100-progress_percentage))/60),mod(t_ave*(100-progress_percentage),60));
            end
        end
        function [y] = countCluster(obj, row, col, oldCounter)
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% Searches the cells which are adjacent (north, east, west and south)
% to the (row,col) cell of obj.A, provided that they are
% within bounds and haven't previously been visited (notice
% that "isSafe" checks for both conditions). Also, update the visited 
% matrix of obj. Returns the number of elements above threshold for
% the cluster.
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
            arguments  
                obj 
                row (1,1) 
                col (1,1)  
                oldCounter (1,1)
            end
            counter = oldCounter;
            adjRow = [-1, 0, 1, 0]; % These arrays are used to find the north, east, 
            adjCol = [0, 1, 0, -1]; % south and west adjacent cells to the current cell.
            threshold = obj.thresholdVal;
            
            if isVisited(obj, row, col) == false % If the cell hasn't been visited...
                visit(obj, row, col); % Set the cell to visited     
            end
            if obj.A(row,col)>= threshold
                counter = counter + 1;
                obj.clusterMatrix(row,col)=1; % Set the element as part of a cluster
                for i=1:4 % Repeat this loop 4 times (4 adjacent cells)
                    if isSafe(obj, row+adjRow(i), col+adjCol(i)) == false
                    end
                    if isSafe(obj, row+adjRow(i), col+adjCol(i)) % If the adjacents cell is safe (within bounds and unvisited)
                        visit(obj, row+adjRow(i), col+adjCol(i)); % Visit the adjacent cells after checking if they"re safe
                        if obj.A(row+adjRow(i),col+adjCol(i))>= threshold % If above threshold, increment counter
                            counter = countCluster(obj, row+adjRow(i), col+adjCol(i), counter); % Updates the counter after recursing
                            countCluster(obj, row+adjRow(i), col+adjCol(i), counter); % Recurse the function
                        end
                    end
                end     
            end
            y = counter;
        end
        function [clusters] = searchClusters(obj,maxSize)
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% Searches the obj.A matrix for clusters. maxSize
% specifies the maximum dimensions of the clusters to be
% searched. This method returns an array "clusters", where
% index "k" equals the number of times that a cluster of dimension 
% k has been found within the matrix. This method does NOT
% filter for invalid clusters.
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
            arguments
                obj       
                maxSize 
            end
            obj.visMatrix = zeros(obj.numCols,obj.numRows);
            obj.clusterMatrix = zeros(obj.numCols,obj.numRows); % These lines are used to reset the matrices
            clusters = zeros(1,maxSize); % This is the array that keeps track of the clusters that are found.
            fprintf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
            fprintf("Initiating cluster counting\n");
            fprintf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
            fprintf("Progress: 0 %% // Elapsed: 0 min, 0.000000 s // Remaining: unknown.\n");
            tic
                for row = (1:obj.numCols)
                    for col = (1:obj.numRows)
                        numCols_approx = obj.numCols-mod(obj.numCols,100); % This is 2000 if numCols = 2048 pixels
                        if mod(row*100,numCols_approx) == 0 && col == 1 && row <= numCols_approx % This occurs every approx 1% progress increment
                            t = toc; % Keeps track of elapsed time
                            progress_percentage = floor((row/floor(numCols_approx/100)));
                            t_ave = t/progress_percentage; % Average time for 1% progress 
                            fprintf("Progress: %d %% // ",progress_percentage);
                            fprintf("Elapsed: %d min, %f s // Remaining: %d min, %f s.\n", floor(t/60),mod(t,60),floor((t_ave*(100-progress_percentage))/60),mod(t_ave*(100-progress_percentage),60));
                        end
                        if (obj.clusterMatrix(row,col) == 0) && (isVisited(obj,row,col) == 0)  % Only if the element isn"t already part of a 
                                                                                               % previously found cluster and unvisited
                            temp = countCluster(obj,row,col,0);
                            if temp <= maxSize && temp ~= 0
                                clusters(temp)=clusters(temp)+1;
                            end
                        end
                    end
                end
                % Print the results of the search after completion  
                fprintf("\nClusters search complete. Total elapsed time: %d min, %f s.\n", floor(t/60),mod(t,60));
                fprintf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
                fprintf("Cluster statistics:\n");
                fprintf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
                for i=1:maxSize
                    fprintf("Number of cluster with dimension %d: %d\n",i,clusters(i));
                end
            end
        function [y1,y2,y3,y4,y_tot,dim1,dim2,dim3,dim4,dim_tot] = makeHistogram(obj,numClusters)
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% Creates a set of histograms: histogram1,histogram2,histogram3,histogram4
% and total histogram (sum of all histograms)
% of the xray spectrum for the selected object.
% numClusters specifies the maximum dimension of the clusters
% to be considered to construct the array. Note that
% numClusters must be an integer between 1 and 4, and the
% clusters are filtered for bad clusters. NOTE: This functions works well
% for 1 image, but if using multiple images it is better to use
% the function "fastCoordClustersHist" (although you need to
% run "findCoordClusters" before that and pass it as an
% input). This function does not return the coordinates of the
% clusters.
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
            arguments
                obj
                numClusters   
            end
            obj.visMatrix = zeros(obj.numCols,obj.numRows);
            obj.clusterMatrix = zeros(obj.numCols,obj.numRows); % These lines are used to reset the matrices
            
            clust1 = 0; % These will count the total clusters used for generating the histogram
            clust2 = 0;
            clust3 = 0;
            clust3_valid = 0;
            clust4 = 0;
            clust4_valid = 0;

            histogram = zeros(1,2^16); % This array will hold the bins of the total histogram
            % These arrays are preallocated for the partial histograms y_i
            histogram1 = histogram; % Partial histogram of single clusters
            histogram2 = histogram; % Partial histogram of double clusters
            histogram3 = histogram; % Partial histogram of triple clusters
            histogram4 = histogram; % Partial histogram of quadruple clusters  
            switch numClusters
                case 1 % Only consider single clusters to be valid
                    for row = (1:obj.numCols)
                        for col = (1:obj.numRows)
                            if (obj.clusterMatrix(row,col) == 0) && (isVisited(obj,row,col) == 0)  % Only if the element isn"t already part of a 
                                                                                                   % previously found cluster and unvisited
                                dimCluster = countCluster(obj,row,col,0);

                                if dimCluster == 1
                                    histogram1(obj.A(row,col)-obj.noiseVal) = histogram1(obj.A(row,col)-obj.noiseVal) + 1;   % Anytime a single cluster is found
                                                                                                                             % Noise is also subtracted here.
                                    clust1 = clust1 + 1; % Count the single clusters for later statistics
                                end
                            end
                        end
                    end
                    y1=histogram1;
                    y2=histogram2;
                    y3=histogram3;
                    y4=histogram4;
                    y_tot = histogram1+histogram2+histogram3+histogram4;
   
                case 2
                    for row = (1:obj.numCols)
                        for col = (1:obj.numRows)
                            if (obj.clusterMatrix(row,col) == 0) && (isVisited(obj,row,col) == 0)  % Only if the element isn"t already part of a 
                                                                                                   % previously found cluster and unvisited
                                dimCluster = countCluster(obj,row,col,0);

                                if dimCluster == 1 || dimCluster == 2
                                    if dimCluster == 1
                                        histogram1(obj.A(row,col)-obj.noiseVal) = histogram1(obj.A(row,col)-obj.noiseVal) + 1;   % Anytime a single cluster is found
                                                                                                                               % increment bin corresponding to the
                                                                                                                               % value of the counts in the location
                                                                                                                               % of the cluster (single photon event).
                                                                                                                               % Noise is also subtracted here.
                                        clust1 = clust1 + 1; % Count the single clusters for later statistics
                                    end                                                                                    
                                    
                                    if dimCluster == 2
                                        pixel_count1 = obj.A(row,col)-obj.noiseVal; % This is only part of the total count, we must add both pixels in the cluster
                                        % The ifs below take into account
                                        % that if a double cluster is
                                        % found, the cell above threshold to be considered is either below or to the right of the
                                        % current cell.
                                        clust2 = clust2 + 1; % Count the double clusters for later statistics
                                        if isInBounds(obj,row+1,col)
                                            if obj.clusterMatrix(row+1,col) == 1  % Below and within bounds
                                                pixel_count2 = obj.A(row+1,col)-obj.noiseVal;
                                            end
                                        end
                                        if isInBounds(obj,row,col+1)
                                            if obj.clusterMatrix(row,col+1) == 1 % To the right and within bounds
                                                pixel_count2 = obj.A(row,col+1)-obj.noiseVal;

                                            end
                                        end
                                        pixel_count12 = pixel_count1 + pixel_count2; % When double cluster is found, add the values up and consider as single event.
                                        histogram2(pixel_count12) = histogram2(pixel_count12) + 1;   % Anytime a single cluster is found
                                                                                       % increment bin corresponding to the
                                                                                       % value of the counts in the location
                                                                                       % of the cluster (single photon event).
                                                                                       % Noise is also subtracted (before)
                                    end                                                                                       
                                end
                            end
                        end
                    end
                    y1=histogram1;
                    y2=histogram2;
                    y3=histogram3;
                    y4=histogram4;
                    y_tot = histogram1+histogram2+histogram3+histogram4;
                    
                case 3
                    for row = (1:obj.numCols)
                        for col = (1:obj.numRows)
                            if (obj.clusterMatrix(row,col) == 0) && (isVisited(obj,row,col) == 0)  % Only if the element isn"t already part of a 
                                                                                                   % previously found cluster and unvisited
                                dimCluster = countCluster(obj,row,col,0);

                                if dimCluster == 1 || dimCluster == 2 || dimCluster == 3
                                    if dimCluster == 1
                                        histogram1(obj.A(row,col)-obj.noiseVal) = histogram1(obj.A(row,col)-obj.noiseVal) + 1;   % Anytime a single cluster is found
                                                                                                                               % increment bin corresponding to the
                                                                                                                               % value of the counts in the location
                                                                                                                               % of the cluster (single photon event).
                                                                                                                               % Noise is also subtracted here.
                                        clust1 = clust1 + 1; % Count the single clusters for later statistics
                                    end                                                                                    
                                    
                                    if dimCluster == 2
                                        clust2 = clust2 + 1; % Count the double clusters for later statistics
                                        pixel_count1 = obj.A(row,col)-obj.noiseVal; % This is only part of the total count, we must add both pixels in the cluster
                                        % The ifs below take into account
                                        % that if a double cluster is
                                        % found, the cell above threshold to be considered is either below or to the right of the
                                        % current cell.
                                        
                                        if isInBounds(obj,row+1,col)
                                            if obj.clusterMatrix(row+1,col) == 1  % Below and within bounds
                                                pixel_count2 = obj.A(row+1,col)-obj.noiseVal;
                                            end
                                        end
                                        if isInBounds(obj,row,col+1)
                                            if obj.clusterMatrix(row,col+1) == 1 % To the right and within bounds
                                                pixel_count2 = obj.A(row,col+1)-obj.noiseVal;

                                            end
                                        end
                                        pixel_count12 = pixel_count1 + pixel_count2; % When double cluster is found, add the values up and consider as single event.
                                        histogram2(pixel_count12) = histogram2(pixel_count12) + 1;   % Anytime a single cluster is found
                                                                                       % increment bin corresponding to the
                                                                                       % value of the counts in the location
                                                                                       % of the cluster (single photon event).
                                                                                       % Noise is also subtracted (before)
                                    end 
                                    if dimCluster == 3
                                         
                                        clust3 = clust3 + 1; % Count the triple clusters for later statistics
                                        pixel_count_temp = 0;
                                        % Once a triple cluster is
                                        % found, we check that it is
                                        % not in a line (which is
                                        % invalid). Then we check the
                                        % cells in this way:
                                        %         ..OOOOO..
                                        %         ..OOCXO..
                                        %         ..OXXXO..
                                        %         ..OOOOO..
                                        % Where C is current cell, O
                                        % are unchecked cells, X are
                                        % checked. 
                                        adjRow = [0, 0, 1, 1, 1]; % These arrays are used to find the cells described above
                                        adjCol = [0, 1, 1, 0, -1];
                                        k = 0; % This is used to check that we aren"t counting some other cluster
                                        for i=1:5    
                                            if isInBounds(obj,row+adjRow(i),col+adjCol(i)) && obj.clusterMatrix(row+adjRow(i),col+adjCol(i)) == 1
                                                pixel_count_temp = pixel_count_temp+obj.A(row+adjRow(i),col+adjCol(i))-obj.noiseVal;
                                                k = k + 1;
                                            end
                                        end
                                        if isInBounds(obj,row,col+2)==true
                                            if k == 3 && (obj.clusterMatrix(row,col+2) ~= 1)
                                                clust3_valid = clust3_valid + 1; % Count the valid triple clusters for later statistics
                                                pixel_count123 = pixel_count_temp; % The triple cluster is L shaped (valid)
                                                histogram3(pixel_count123) = histogram3(pixel_count123) + 1;
                                            end
                                        elseif isInBounds(obj,row,col+2)==false
                                            if k == 3
                                                clust3_valid = clust3_valid + 1; % Count the valid triple clusters for later statistics
                                                pixel_count123 = pixel_count_temp; % The triple cluster is L shaped (valid)
                                                histogram3(pixel_count123) = histogram3(pixel_count123) + 1;
                                            end     
                                        end  
                                    end
                                end
                            end
                        end
                    end
                    fprintf("\nSingle, double and triple clusters histogram complete.\n"); 
                    fprintf("There are %d triple clusters, of which %d are valid (%f%%).\n",clust3,clust3_valid,(clust3_valid/clust3)*100);
                    y1=histogram1;
                    y2=histogram2;
                    y3=histogram3;
                    y4=histogram4;
                    y_tot = histogram1+histogram2+histogram3+histogram4;
                    
                case 4
                    for row = (1:obj.numCols)
                        for col = (1:obj.numRows)
                            if (obj.clusterMatrix(row,col) == 0) && (isVisited(obj,row,col) == 0)  % Only if the element isn"t already part of a 
                                                                                                   % previously found cluster and unvisited
                                dimCluster = countCluster(obj,row,col,0);
                                if dimCluster == 1 || dimCluster == 2 || dimCluster == 3 || dimCluster == 4
                                    if dimCluster == 1
                                        histogram1(obj.A(row,col)-obj.noiseVal) = histogram1(obj.A(row,col)-obj.noiseVal) + 1;   % Anytime a single cluster is found
                                                                                                                               % increment bin corresponding to the
                                                                                                                               % value of the counts in the location
                                                                                                                               % of the cluster (single photon event).
                                                                                                                               % Noise is also subtracted here.
                                        clust1 = clust1 + 1; % Count the single clusters for later statistics
                                    end                                                                                    
                                    
                                    if dimCluster == 2
                                        clust2 = clust2 + 1; % Count the double clusters for later statistics
                                        pixel_count1 = obj.A(row,col)-obj.noiseVal; % This is only part of the total count, we must add both pixels in the cluster
                                        % The ifs below take into account
                                        % that if a double cluster is
                                        % found, the cell above threshold to be considered is either below or to the right of the
                                        % current cell.
                                        
                                        if isInBounds(obj,row+1,col)
                                            if obj.clusterMatrix(row+1,col) == 1  % Below and within bounds
                                                pixel_count2 = obj.A(row+1,col)-obj.noiseVal;
                                            end
                                        end
                                        if isInBounds(obj,row,col+1)
                                            if obj.clusterMatrix(row,col+1) == 1 % To the right and within bounds
                                                pixel_count2 = obj.A(row,col+1)-obj.noiseVal;

                                            end
                                        end
                                        pixel_count12 = pixel_count1 + pixel_count2; % When double cluster is found, add the values up and consider as single event.
                                        histogram2(pixel_count12) = histogram2(pixel_count12) + 1;   % Anytime a single cluster is found
                                                                                       % increment bin corresponding to the
                                                                                       % value of the counts in the location
                                                                                       % of the cluster (single photon event).
                                                                                       % Noise is also subtracted (before)
                                    end 
                                    if dimCluster == 3
                                         
                                        clust3 = clust3 + 1; % Count the triple clusters for later statistics
                                        %pixel_count_temp = obj.A(row,col)-obj.noiseVal;
                                        pixel_count_temp = 0;
                                        % Once a triple cluster is
                                        % found, we check that it is
                                        % not in a line (which is
                                        % invalid). Then we check the
                                        % cells in this way:
                                        %         ..OOOOO..
                                        %         ..OOCXO..
                                        %         ..OXXXO..
                                        %         ..OOOOO..
                                        % Where C is current cell, O
                                        % are unchecked cells, X are
                                        % checked.
                                        adjRow = [0, 0, 1, 1, 1]; % These arrays are used to find the cells described above
                                        adjCol = [0, 1, 1, 0, -1];
                                        k = 0; % This is used to check that we aren"t counting some other cluster
                                        for i=1:5
                                            if isInBounds(obj,row+adjRow(i),col+adjCol(i)) && obj.clusterMatrix(row+adjRow(i),col+adjCol(i)) == 1
                                                pixel_count_temp = pixel_count_temp+obj.A(row+adjRow(i),col+adjCol(i))-obj.noiseVal;
                                                k = k + 1;
                                            end
                                        end
                                        if isInBounds(obj,row,col+2)==true
                                            if k == 3 && (obj.clusterMatrix(row,col+2) ~= 1)
                                                clust3_valid = clust3_valid + 1; % Count the valid triple clusters for later statistics
                                                pixel_count123 = pixel_count_temp; % The triple cluster is L shaped (valid)
                                                histogram3(pixel_count123) = histogram3(pixel_count123) + 1;
                                            end
                                        elseif isInBounds(obj,row,col+2)==false
                                            if k == 3
                                                clust3_valid = clust3_valid + 1; % Count the valid triple clusters for later statistics
                                                pixel_count123 = pixel_count_temp; % The triple cluster is L shaped (valid)
                                                histogram3(pixel_count123) = histogram3(pixel_count123) + 1;
                                            end                                       
                                        end
                                    end
                                    if dimCluster == 4
                                        clust4 = clust4 + 1; % Count the quadruple clusters for later statistics
                                        pixel_count_temp = 0;
                                        % Once a quadruple cluster is
                                        % found, we check that it is
                                        % a square, which is the only 
                                        % valid shape.
                                        % Then we check the
                                        % cells in this way:
                                        %         ..OOOOO..
                                        %         ..OOCXO..
                                        %         ..OOXXO..
                                        %         ..OOOOO..
                                        % Where C is current cell, O
                                        % are unchecked cells, X are
                                        % checked.
                                        adjRow = [0, 0, 1, 1]; % These arrays are used to find the cells described above
                                        adjCol = [0, 1, 1, 0];
                                        k = 0; % This is used to check that we aren"t counting some other cluster
                                        for i=1:4 
                                            if isInBounds(obj,row+adjRow(i),col+adjCol(i)) && obj.clusterMatrix(row+adjRow(i),col+adjCol(i)) == 1
                                                pixel_count_temp = pixel_count_temp+obj.A(row+adjRow(i),col+adjCol(i))-obj.noiseVal;
                                                k = k + 1;
                                            end
                                        end                                
                                        if k == 4 
                                            clust4_valid = clust4_valid + 1; % Count the valid quadruple clusters for later statistics
                                            pixel_count1234 = pixel_count_temp; % The quadruple cluster is square shaped (valid)
                                            histogram4(pixel_count1234) = histogram4(pixel_count1234) + 1;
                                        end  
                                    end
                                end
                            end
                        end
                    end
                    y1=histogram1;
                    y2=histogram2;
                    y3=histogram3;
                    y4=histogram4;
                    y_tot = histogram1+histogram2+histogram3+histogram4;      
            end
            dim1=clust1;
            dim2=clust2;
            dim3=clust3_valid;
            dim4=clust4_valid;
            dim_tot=clust1+clust2+clust3_valid+clust4_valid;
            fprintf("\n---------------------------------\nHistogram of file %c%s%c complete!\nIt contains data from: \n%d single clusters\n%d double clusters\n%d valid triple clusters (+%d invalid, %.2f%% of total)\n%d valid quadruple clusters (+%d invalid, %.2f%% of total)\nA total of: %d photons\n---------------------------------\n",'"',obj.fileName,'"',clust1,clust2,clust3_valid,(clust3-clust3_valid),(clust3-clust3_valid)/clust3*100,clust4_valid,(clust4-clust4_valid),(clust4-clust4_valid)/clust4*100,clust1+clust2+clust3+clust4);
        end
        function [coordClusters] = findCoordClusters(obj)
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% Creates the array "coordClusters", which contains the coordinates of all
% the valid clusters that were identified in the
% obj.A matrix. coordClusters is of this form:
% Size | row | col | Counts | 
%   .  |  .  |  .  |    .   |
%   .  |  .  |  .  |    .   |
%   .  |  .  |  .  |    .   |
% Where total number of rows (minus 1, because of an empty row 
% of zeros) identifies the total number of
% valid clusters, and row, col are the coordinates of the
% pixel of the cluster with greatest value. This function does NOT
% return the histogram of the xray matrix object, but only the
% coordinate of the cluster. Therefore, it has to be coupled
% with the function coordClustersHist for this purpose. The
% function also modifies the field obj.coordCluster at the end,
% so that the cluster coordinates are marked onto the object.
% Note that obj needs to be a SINGLE OBJECT, not an array of
% objs.
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
            arguments
                obj  
            end
            obj.visMatrix = zeros(obj.numRows,obj.numCols);
            obj.clusterMatrix = zeros(obj.numRows,obj.numCols); % These lines are used to reset the matrices
            
            clust1 = 0; % These will count the total clusters used for generating the histogram
            clust2 = 0;
            clust3 = 0;
            clust3_valid = 0;
            clust4 = 0;
            clust4_valid = 0;
            coordClusters = zeros(1,4); % We initialize the coordClusters matrix to this value.
                                        % Afterwards, we will append new
                                        % rows every time a cluster is
                                        % found. If no cluster is found we
                                        % will just get the first row of 
                                        % [0 0 0 0]...  
            fprintf("\nSearching cluster coordinates of file %c%s%c.\n",'"',obj.fileName,'"')
            for row = (1:obj.numCols)
                for col = (1:obj.numRows)
                    if (obj.clusterMatrix(row,col) == 0) && (isVisited(obj,row,col) == 0)  % Only if the element isn"t already part of a 
                                                                                           % previously found cluster and unvisited
                        dimCluster = countCluster(obj,row,col,0);

                        if dimCluster == 1 || dimCluster == 2 || dimCluster == 3 || dimCluster == 4
                            if dimCluster == 1
                                clust1 = clust1 + 1; % Count the single clusters for later statistics
                                coordClusters(end+1,:)=[1,row,col,obj.A(row,col)-obj.noiseVal];
                            end                                                                                    
                            if dimCluster == 2
                                clust2 = clust2 + 1; % Count the double clusters for later statistics
                                pixel_count1 = obj.A(row,col)-obj.noiseVal; % This is only part of the total count, we must add both pixels in the cluster
                                % The ifs below take into account
                                % that if a double cluster is
                                % found, the cell above threshold to be considered is either below or to the right of the
                                % current cell.
                                if isInBounds(obj,row+1,col)
                                    if obj.clusterMatrix(row+1,col) == 1  % Below and within bounds
                                        pixel_count2 = obj.A(row+1,col)-obj.noiseVal;
                                        % --------------- delete this part
                                        % if it doesn"t work
                                        if max(pixel_count1,pixel_count2)==pixel_count2
                                            row_max = row+1;
                                            col_max = col;
                                        else
                                            row_max = row;
                                            col_max = col;
                                        end
                                        % ---------------
                                    end
                                end
                                if isInBounds(obj,row,col+1)
                                    if obj.clusterMatrix(row,col+1) == 1 % To the right and within bounds
                                        pixel_count2 = obj.A(row,col+1)-obj.noiseVal;
                                        % --------------- delete this part
                                        % if it doesn"t work
                                        if max(pixel_count1,pixel_count2)==pixel_count2
                                            row_max = row;
                                            col_max = col+1;
                                        else
                                            row_max = row;
                                            col_max = col;
                                        end
                                        % ---------------
                                    end
                                end
                                pixel_count12 = pixel_count1 + pixel_count2; % When double cluster is found, add the values up and consider as single event.
                                coordClusters(end+1,:)=[2,row_max,col_max,pixel_count12]; % This actually stores the coordinate of the largest element
                                                                                          % of the cluster, instead of the one on the top left corner.
                            end 
                            if dimCluster == 3
                                clust3 = clust3 + 1; % Count the triple clusters for later statistics
                                pixel_count_temp = 0;
                                % Once a triple cluster is
                                % found, we check that it is
                                % not in a line (which is
                                % invalid). Then we check the
                                % cells in this way:
                                %         ..OOOOO..
                                %         ..OOCXO..
                                %         ..OXXXO..
                                %         ..OOOOO..
                                % Where C is current cell, O
                                % are unchecked cells, X are
                                % checked.
                                adjRow = [0, 0, 1, 1, 1]; % These arrays are used to find the cells described above
                                adjCol = [0, 1, 1, 0, -1];
                                k = 0; % This is used to check that we aren"t counting some other cluster
                                row_max = row; % This allows to try to find what pixel has the highest value, initially we set the max value coordinates to those the current pixel
                                col_max = col;
                                for i=1:5
                                    if isInBounds(obj,row+adjRow(i),col+adjCol(i)) && obj.clusterMatrix(row+adjRow(i),col+adjCol(i)) == 1
                                        pixel_count_temp = pixel_count_temp+obj.A(row+adjRow(i),col+adjCol(i))-obj.noiseVal;
                                        if i > 1 && isInBounds(obj,row+adjRow(i-1),col+adjCol(i-1)) % Ensures that we can look at the previous element without going off the indices
                                            if obj.A(row+adjRow(i),col+adjCol(i))>obj.A(row+adjRow(i-1),col+adjCol(i-1)) 
                                                % As we scan through the 5
                                                % elements, if the next element
                                                % has a pixel value greater
                                                % than the one before, then we
                                                % set the coordinates row_max
                                                % and col_max to this new
                                                % value. This is basically
                                                % similar to a bubble sort but
                                                % only for the max element with
                                                % a linear search.
                                                row_max = row+adjRow(i);
                                                col_max = col+adjCol(i);
                                            end
                                        end   
                                        k = k + 1;
                                    end
                                end
                                if isInBounds(obj,row,col+2)==true
                                    if k == 3 && (obj.clusterMatrix(row,col+2) ~= 1)
                                        clust3_valid = clust3_valid + 1; % Count the valid triple clusters for later statistics
                                        pixel_count123 = pixel_count_temp; % The triple cluster is L shaped (valid)
                                        %coordClusters(end+1,:)=[3,row,col,pixel_count123];
                                        coordClusters(end+1,:)=[3,row_max,col_max,pixel_count123];
                                    end
                                elseif isInBounds(obj,row,col+2)==false
                                    if k == 3
                                        clust3_valid = clust3_valid + 1; % Count the valid triple clusters for later statistics
                                        pixel_count123 = pixel_count_temp; % The triple cluster is L shaped (valid)
                                        %coordClusters(end+1,:)=[3,row,col,pixel_count123];
                                        coordClusters(end+1,:)=[3,row_max,col_max,pixel_count123];
                                    end                                       
                                end
                            end
                            if dimCluster == 4
                                clust4 = clust4 + 1; % Count the quadruple clusters for later statistics
                                pixel_count_temp = 0;
                                % Once a quadruple cluster is
                                % found, we check that it is
                                % a square, which is the only 
                                % valid shape.
                                % Then we check the
                                % cells in this way:
                                %         ..OOOOO..
                                %         ..OOCXO..
                                %         ..OOXXO..
                                %         ..OOOOO..
                                % Where C is current cell, O
                                % are unchecked cells, X are
                                % checked.
                                adjRow = [0, 0, 1, 1]; % These arrays are used to find the cells described above
                                adjCol = [0, 1, 1, 0];
                                k = 0; % This is used to check that we aren"t counting some other cluster
                                for i=1:4
                                    if isInBounds(obj,row+adjRow(i),col+adjCol(i)) && obj.clusterMatrix(row+adjRow(i),col+adjCol(i)) == 1
                                        pixel_count_temp = pixel_count_temp+obj.A(row+adjRow(i),col+adjCol(i))-obj.noiseVal;
                                        if i > 1 && isInBounds(obj,row+adjRow(i-1),col+adjCol(i-1)) % Ensures that we can look at the previous element without going off the indices
                                            if obj.A(row+adjRow(i),col+adjCol(i))>obj.A(row+adjRow(i-1),col+adjCol(i-1)) 
                                                % As we scan through the 4
                                                % elements, if the next element
                                                % has a pixel value greater
                                                % than the one before, then we
                                                % set the coordinates row_max
                                                % and col_max to this new
                                                % value. This is basically
                                                % similar to a bubble sort but
                                                % only for the max element with
                                                % a linear search.
                                                row_max = row+adjRow(i);
                                                col_max = col+adjCol(i);
                                            end
                                        end
                                        k = k + 1;
                                    end
                                end                                
                                if k == 4 
                                    clust4_valid = clust4_valid + 1; % Count the valid quadruple clusters for later statistics
                                    pixel_count1234 = pixel_count_temp; % The quadruple cluster is square shaped (valid)
                                    coordClusters(end+1,:)=[4,row_max,col_max,pixel_count1234];
                                end  
                            end
                        end
                    end
                end
            end
            if size(coordClusters,1)==1 % No clusters have been counted, return a [0 0 0 0]
                coordClusters = [0 0 0 0]; % This is probably redundant, but just to make sure...
            elseif size(coordClusters,1)>1 % Some clusters have been counted, therefore eliminate the empty first row
                coordClusters(1,:)=[]; 
            end
            obj.coordClusters = coordClusters;
            fprintf("\n---------------------------------\nCluster coordinates of file %c%s%c found!\nData from: \n%d single clusters\n%d double clusters\n%d triple clusters\n%d quadruple clusters\nA total of: %d photons\n---------------------------------\n",'"',obj.fileName,'"',clust1,clust2,clust3_valid,clust4_valid,clust1+clust2+clust3_valid+clust4_valid);  
        end
        function [sum_hist, x_hist] = fastCoordClustersHist(obj,xmin,xmax,bin_size,rect_window)
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% Use this function to create a histogram of a multiple matrix 
% objects in an array (could also be just one), using a 
% parallel "parfor" loop. The input for this function is a 1xN
% vector of xray matrix objects, a lower limit for energy to be
% considered (xmin), a higher limit for energy to be considered
% (xmax), the bin size of the histogram in eV (bin_size) and a
% Mx4 matrix "rect_windows". This matrix provides the
% coordinates of rectangles that will be used to selectively
% calculate the histogram. The matrix is of this form:
%
%           1_r1  |  1_c1  |  1_r2  |  1_c2  |
%           2_r1  |  2_c1  |  2_r2  |  2_c2  |
%             .   |    .   |    .   |    .   |
%             .   |    .   |    .   |    .   |
%           N_r1  |  N_c1  |  N_r2  |  N_c2  |

% Where i_r1 is row and i_c1 is the column of the upper-left
% corner of a rectangle "i" that's within the dimensions of the
% xray object matrix, and i_r2 is row and i_c2 is the column of 
% the bottom-right corner of that same rectangle "i". This
% means that the entire matrix stores the coordinates of M
% different rectangular windows to examine. Make sure that each
% rectangle has coordinates that do not exceed the dimensions
% of the matrix upon which they should be applied!
% The function returns the array sum_hist, which is a 1xM
% row vector whose ith element is the sum of the histograms of
% N = size(obj,2) images calculated over the ith rectangular
% window. 
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
            arguments
                obj
                xmin
                xmax
                bin_size
                rect_window
            end
            edges = xmin:bin_size:xmax;
            x_hist = 0.5*(edges(1:end-1)+edges(2:end)); % We want to plot the amount of counts for each bin in the center of the bin, which we can obtain by taking this average of the edges...
            sum_hist = zeros(size(rect_window,1),size(x_hist,2));
            merged_cluster = [0 0 0 0];
            N = size(obj,2); % Number of images
            M = size(rect_window,1); % Number of windows
            tic;
            parfor i = 1:N % Merge all the cluster coordinates into one big matrix for all the N images
                merged_cluster = [merged_cluster;obj(i).coordClusters];
            end
            toc
            merged_cluster(1,:)=[]; % We initialized merged_cluster to [0 0 0 0] so we now delete that row
            C = merged_cluster; % Just done for ease of readability of the following...
            number_Rows = obj(1).numRows;
            number_Cols = obj(1).numCols;
            for j = 1:M % For as many rows as rect_window has (the number of windows we have, M)
                fprintf("Working on window %d of %d\n",j,M)
                hist_temp = zeros(1,size(x_hist,2));
                r1 = rect_window(j,1);
                c1 = rect_window(j,2);
                r2 = rect_window(j,3);
                c2 = rect_window(j,4);
                coordClustersCrop = C((C(:,2)>=r1)&(C(:,2)<=r2)&(C(:,3)>=c1)&(C(:,3)<=c2),:);
                [hist_curr,~]=histcounts(coordClustersCrop(:,4)*11.7,edges); % Here we are creating the total histogram in eV already with the conversion that 1 count = 11.7 eV
                hist_temp = hist_temp + hist_curr;
                sum_hist(j,:) = hist_temp;
            end
        end
    end
end