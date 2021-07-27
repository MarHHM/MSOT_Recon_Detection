function [base_spectra] = loadSpectra_m(spectra_path, wavelengths, INCLUDE_BLOOD, NORMALIZE_BASE_SPECTRA)
% - VIP: this function doesn't load any endogenous chromohpores spectra other than oxy & deoxy blood (even if the spectrum file of this chromohpore exists in the spectra_path)
% - spectra -> n*l matrix (n: number of agents + oxy + deoxy - l: len(wavelengths))
% - generated spectra are spline-interpolation of the full spectra (coming from spectrometer) for the acquisition wavelengths

keep_path = pwd;

cd(spectra_path);
if (nargin<2)
  load('wavelengths.mat');
end

n_wavelengths = length(wavelengths);
generaldirectory1 = dir ('*agent_*');
n_spectra = size(generaldirectory1, 1);
if INCLUDE_BLOOD
  generaldirectory2 = dir ('*hb_spectra*');
  if ~size(generaldirectory2, 1)
    error("No blood spectra file found (hb_spectra.mat)!");
  end
  n_spectra = n_spectra + 2;             % 2 --> deoxy & oxy spectra
end

base_spectra = zeros(n_spectra, n_wavelengths);
base_spectra_names = strings(n_spectra, 1);

k=1;
%load all the agents (row-by-row) first in the folder
for i=1:size(generaldirectory1,1)
  load(generaldirectory1(i).name);
  agent_name = generaldirectory1(i).name(7:end-4);
  base_spectra_names(k) = agent_name;
  var = eval(agent_name);
  base_spectra(k,:) = spline(var(:,1),var(:,2),wavelengths);
  k=k+1;
end

%then load deoxy then oxy blood spectra
if INCLUDE_BLOOD
  load(generaldirectory2(1).name);
  
  %Calibrate, at 800 nm both should be equal to 1.
  idx = find(abs(hg_spectra(:,1) - 800) == min(abs(hg_spectra(:,1) - 800)));
  hg_spectra(:,2) = hg_spectra(:,2)./hg_spectra(idx,2);
  hg_spectra(:,3) = hg_spectra(:,3)./hg_spectra(idx,3);
  
  %Resample
  base_spectra(k,:) = spline(hg_spectra(:,1),hg_spectra(:,3),wavelengths);
  base_spectra_names(k) = "deoxy blood";
  k = k + 1;
  base_spectra(k,:) = spline(hg_spectra(:,1),hg_spectra(:,2),wavelengths);
  base_spectra_names(k) = "oxy blood";
end
disp("Base spectra for unmixing are (in order): ")
for i=1:size(base_spectra_names)
  fprintf("\t" + base_spectra_names(i) + "\n")
end

if NORMALIZE_BASE_SPECTRA
  for i=1:size(base_spectra,1)
    base_spectra(i,:) = base_spectra(i,:)./norm(base_spectra(i,:));
  end
end

cd(keep_path);

end