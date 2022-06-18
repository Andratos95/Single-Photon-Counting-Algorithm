function [hist_tot,x_hist,dim] = coordClustersAnalyzer (coordClusters,xmin,xmax,bin_size)
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% Returns the histogram of the coordClusters list in the range between
% [xmin,xmax] (in eV) and at the bin size of bin_size (in eV also). It also
% returns the amount of VALID clusters of dimension 1,2,3 and 4, along with 
% the total number of clusters in a row vector "dim" where:
% dim = [dim1 dim2 dim3 dim4 dim_tot]
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %

    arguments
        coordClusters % A list of coordinates of all the clusters
        xmin % The initial range of energies to be used in eV (often 2 keV or so...)
        xmax % The final range of energies to be used in eV (often 25 keV or so...)
        bin_size % The size of the bins for the histogram in eV (should be equal to RMS of noise in eV)
    end
    
    dim1 = size(coordClusters(coordClusters(:,1)==1,:),1); % Logical indexing is used to select only the rows of coordClusters for which the first element is 1 (2,3,4 for subsequent lines)
    dim2 = size(coordClusters(coordClusters(:,1)==2,:),1);
    dim3 = size(coordClusters(coordClusters(:,1)==3,:),1);
    dim4 = size(coordClusters(coordClusters(:,1)==4,:),1);
    dim_tot = dim1+dim2+dim3+dim4;
    dim = [dim1 dim2 dim3 dim4 dim_tot];
    
    [hist_tot,edges]=histcounts(coordClusters(:,4)*11.7,xmin:bin_size:xmax); % Here we are creating the total histogram in eV already with the conversion that 1 count = 11.7 eV
    x_hist=0.5*(edges(1:end-1)+edges(2:end)); % We want to plot the amount of counts for each bin in the center of the bin, which we can obtain by taking this average of the edges...


    
end