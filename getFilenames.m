function y = getFilenames(date, initial_filename, num_images)
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% Gets the filenames of the images from a certain date,
% starting from the initial_filename image and up to the
% successive num_images (can be passed to the function 
% sum_xray_histograms). Returns a vector 'filenames' which
% contains num_images strings corresponding to the file names
% of the images we want to analyze.
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
            arguments
                date 
                initial_filename
                num_images
            end
  
            dir_name = strcat('./x-rays/',date,'/'); % Goes in the x-rays folder on the correct date
            pat = digitsPattern;
            filename_numbers = extract(initial_filename,pat)';
            batch_number_str = filename_numbers(1); % This is the batch number as string (it has zero padding, like '0002')
            shot_number = str2num(filename_numbers(2));  % This is the shot number as an integer (warning: no zero padding!!! Example: '0112' turns into 112)

            for i = 1:num_images

                % First, we check if we need to add zero padding to the shot number
                if ceil(log10(max(1,abs(shot_number)+1)))<4 % Counts the number of digits of shot_number to add zero padding to the string
                    num_zeros_padding = 4 - ceil(log10(max(1,abs(shot_number)+1))); % this is how many zeros we need to add at the beginning, for padding
                    shot_number_str = num2str(shot_number);
                    for k = 1:num_zeros_padding
                        shot_number_str = strcat('0', shot_number_str); % Pads the necessary zeros at the beginning
                    end   
                end
                % Second, we save the file name into a string array 'filenames'
                % which contains all the filenames of the images we need to analyze

                filenames(i) = strcat(dir_name,batch_number_str,'_',shot_number_str,'_X-ray.tif');
                %sprintf('%s\n',filenames(i))
                shot_number = shot_number + 1;
            end
            y = filenames;
        end