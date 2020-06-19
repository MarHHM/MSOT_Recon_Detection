function [datainfo, sigMat, sigMat_pathName] = loadSigMat_iThera(scan_path)

% [FileNamesigMat, sigMat_pathName, ~] = uigetfile({'*.msot'}, 'mytitle', scan_path);
% datainfo = loadMSOT([sigMat_pathName '\' FileNamesigMat]);

sigMat_pathName = scan_path;
scan_path_parts = strsplit(scan_path, '\');
datainfo = loadMSOT(strcat(scan_path, '\', scan_path_parts(end), ".msot"));
sigMat = loadMSOTSignals(datainfo);
