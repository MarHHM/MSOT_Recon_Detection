function [spectra] = LoadSpectra(data_path, wavelengths) 

keep_path = pwd;

cd(data_path);
if (nargin<2)
    load('wavelengths.mat');
end

 nb_wavelengths = length(wavelengths);
 generaldirectory1 = dir ('*agent_*');
 generaldirectory2 = dir ('*hb_spectra*');
 
 if size(generaldirectory2,1)
    Nb_spectra = size(generaldirectory1,1) + 2;
 else
    Nb_spectra = size(generaldirectory1,1);
 end
    
 spectra = ones(Nb_spectra,nb_wavelengths);
 k=1;
 %load all the agents in the folder
 for i=1:size(generaldirectory1,1)
     load(generaldirectory1(i).name);
     var = eval(generaldirectory1(i).name(7:end-4));
     spectra(k,:) = spline(var(:,1),var(:,2),wavelengths);
     k=k+1;
 end
  
 %Load oxy - deoxy
 if size(generaldirectory2,1)
        load(generaldirectory2(1).name);
        
        %Callibrate, at 800 nm both should be equal to 1.
        idx = find(abs(hg_spectra(:,1) - 800) == min(abs(hg_spectra(:,1) - 800)));
        hg_spectra(:,2) = hg_spectra(:,2)./hg_spectra(idx,2);
        hg_spectra(:,3) = hg_spectra(:,3)./hg_spectra(idx,3);
        
        %Resample
        spectra(k,:) = spline(hg_spectra(:,1),hg_spectra(:,3),wavelengths);
        spectra(k+1,:) = spline(hg_spectra(:,1),hg_spectra(:,2),wavelengths);
        clear hg_spectra;
 end

%  Normalization can only be applied when we are looking for agents and not in SO2 estimation. 
%  for i=1:size(spectra,1)
%     spectra(i,:) = spectra(i,:)./norm(spectra(i,:)); 
%  end
%  
 cd(keep_path);
 
end