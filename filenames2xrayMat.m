function y = filenames2xrayMat (fileNames)
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% Input a vector of filenames and generate a vector of xrayMatrix objects.
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
    for i=1:size(fileNames,2)
        y_temp(i)=xrayMatrix(fileNames(i));
    end
    y = y_temp;
    
end
