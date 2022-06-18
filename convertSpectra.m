function spectra = convertSpectra (hist,x_hist)
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% This function converts the spectra from the matrix "hist" (each row is the
% histogram of different image or part of image) to take into account the
% quantum efficiency QE and the transmission through 3 microns of aluminum,
% 50 microns of kapton tape, 250 microns of beryllium and 
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
arguments
    hist
    x_hist
end
% The quantum efficiency of the BR-DD sensor of the CCD camera

q_eff = table2array(readtable('BR-DD_QE.xlsx','Range','B8:C92')); 
q_eff_interp=interp1(q_eff(:,1),q_eff(:,2),x_hist); % Linear interpolation of the QE on the energies vector
q_eff_interp=q_eff_interp/100; % The QE is given in percentage...
q_ineff_interp=1-q_eff_interp;

% ALUMINUM

tx_Al = table2array(readtable("tx_aluminium_3um.dat"));
tx_Al_interp=interp1(tx_Al(:,1),tx_Al(:,2),x_hist); % Linear interpolation of the aluminum transmission on the energies vector
abs_Al_interp = 1-tx_Al_interp;

% KAPTON

tx_Kapton = table2array(readtable("tx_kapton_50um.dat"));
tx_Kapton_interp=interp1(tx_Kapton(:,1),tx_Kapton(:,2),x_hist); % Linear interpolation of the kapton transmission on the energies vector
abs_Kapton_interp = 1-tx_Kapton_interp;

% BERYLLIUM

tx_Be = table2array(readtable("tx_beryllium_250um.dat"));
tx_Be_interp=interp1(tx_Be(:,1),tx_Be(:,2),x_hist); % Linear interpolation of the kapton transmission on the energies vector
abs_Be_interp = 1-tx_Be_interp;

% CONVERT BY DIVISION

tx_Be_interp=ones(1,size(hist,1));
tx_Kapton_interp=ones(1,size(hist,1));
tx_Al_interp=ones(1,size(hist,1));
for i = 1:size(hist,1)
    spectra_temp(i,:) = hist(1,:)./q_eff_interp./tx_Al_interp./tx_Kapton_interp./tx_Be_interp;
end

spectra = spectra_temp;
 

end