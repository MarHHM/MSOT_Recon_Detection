function [datainfo, Recon_MB] = loadRecon_iThera(scanPath)

[FileNameRecon_MB, ReconMB_pathName, ~] = uigetfile({'*.msot'}, 'mytitle', scanPath);
datainfo = loadMSOT([ReconMB_pathName '\' FileNameRecon_MB]);
Recon_MB = loadMSOTRecon(datainfo);
