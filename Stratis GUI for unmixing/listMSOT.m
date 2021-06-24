function datainfoarr  = listMSOT( cfile )
% LISTMSOT  List .msot META information from Study Folder
%   datainfo = LISTMSOT()      List scans in current directory
%   datainfo = LISTMSOT(dir)   Fist scans in the specified directory
% 
% This function lists the META information on the selected study folder
%
% Required Parameters:
%   <none>
%
% Optional:
%   1: dirname:    path of the study to load.
%
% Return Values:
%   1: datainfoarr: a structure containing the meta information, [] on error 
%
% Example:
%   datainfo = loadMSOT('Scan_1\Scan_1.msot');


%% parameters

xml_filename = cfile;
bsl = strfind(xml_filename,'\');
dirname = xml_filename(1:bsl(end)-1);

try
    patdoc = javaMethod( 'parse', 'com.itheramedical.msotbeans.DataModelMsotProjectDocument$Factory', java.io.File( xml_filename ) ) ;
catch ex
    error( ['Selected XML file could not be loaded. Check well-formedness.\n' ex.message]) ;
end ;
dm = patdoc.getDataModelMsotProject() ;

datainfoarr.DirName = char(dm.getFriendlyName);
datainfoarr.FileName = xml_filename;
datainfoarr.CreationTimeTxt = char(dm.getCreationTime.toString);
datainfoarr.CreationTime = datenum(datainfoarr.CreationTimeTxt,'yyyy-mm-ddTHH:MM:SS.FFF');
datainfoarr.Name = char(dm.getScanNode.getName);
datainfoarr.Comment = char(dm.getScanNode.getComment);
datainfoarr.RealPath = dirname;
datainfoarr.Complete = false;

hw = dm.getHARDWAREDESC;
datainfoarr.HWDesc.TransducerType = char(hw.getTRANSDUCER);
datainfoarr.is3D = strcmp(datainfoarr.HWDesc.TransducerType,'msot3');
datainfoarr.is2D = strcmp(datainfoarr.HWDesc.TransducerType,'msot2');
if (datainfoarr.is3D),
    datainfoarr.TType = '3D';
else
    datainfoarr.TType = '2D';
end

md = dm.getMEASUREMENTDESC;
datainfoarr.USpresent = false;
try
    datainfoarr.USpresent = logical(md.getULTRASOUNDPRESENT);
end


rna = dm.getReconNodes.getDataModelNewReconstructionNodeArray;
datainfoarr.ReconNode = [];
for rn = 1:length(rna)
   rnode = rna(rn);
   datainfoarr.ReconNode(rn).Name = char(rnode.getName);
   datainfoarr.ReconNode(rn).Comment = char(rnode.getComment);
   datainfoarr.ReconNode(rn).GUID = char(rnode.getGUID);
   datainfoarr.ReconNode(rn).Method = char(rnode.getMethod);
   datainfoarr.ReconNode(rn).Resolution = rnode.getResolution;
   datainfoarr.ReconNode(rn).Projections = rnode.getProjections;
   datainfoarr.ReconNode(rn).ROI = double(rnode.getRoi);
end


mna = dm.getMspNodes.getDataModelNewMspNodeArray;
for mn = 1:length(mna)
    mnode = mna(mn);
    datainfoarr.MSPNode(mn).Name = char(mnode.getName);
    datainfoarr.MSPNode(mn).Comment = char(mnode.getComment);
    datainfoarr.MSPNode(mn).Method = char(mnode.getMethod);
    datainfoarr.MSPNode(mn).GUID = char(mnode.getGUID);
    datainfoarr.MSPNode(mn).ReconGUID = char(mnode.getReconGUID);
    % Related Wavelengths
    wla = mnode.getRelatedWavelengths.getWavelengthArray;
    for wl = 1:numel(wla)
        datainfoarr.MSPNode(mn).Wavelengths(wl) = wla(wl).getIntValue;
    end
    % Input Spectra
    datainfoarr.MSPNode(mn).InputSpectra = cell(mnode.getInputSpectra.getStringArray);
    % Recon Info
    rn = 1;
    while ~strcmp(datainfoarr.ReconNode(rn).GUID,datainfoarr.MSPNode(mn).ReconGUID)
        rn = rn + 1;
        if (rn > numel(datainfoarr.ReconNode)) rn = 0; break; end;
    end
    datainfoarr.MSPNode(mn).ReconNodeID = rn;
    if (rn > 0)
        datainfoarr.MSPNode(mn).ReconMethod = datainfoarr.ReconNode(rn).Method;
        datainfoarr.MSPNode(mn).Resolution = datainfoarr.ReconNode(rn).Resolution;
        datainfoarr.MSPNode(mn).Projections = datainfoarr.ReconNode(rn).Projections;
        datainfoarr.MSPNode(mn).ROI = datainfoarr.ReconNode(rn).ROI;
    else
        warning(['Recon ' datainfoarr.MSPNode(mn).ReconGUID ' not found.\n']);
    end
end





%%
% if numel(datainfoarr) > 0,
%     [tmp ind] = sort(arrayfun(@(x) x.CreationTime, datainfoarr));
%     datainfoarr = datainfoarr;
% 
%     for ind = 1:numel(datainfoarr)
%         fprintf('%s\t\t%s\t%s\t%s\n',datainfoarr.DirName,datainfoarr.CreationTimeTxt,datainfoarr.TType,datainfoarr.Name);
%     end
% end

return


