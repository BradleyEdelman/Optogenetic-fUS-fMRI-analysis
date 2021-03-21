function Read_Bruker(inputfolderlist)
% Read Bruker 2dseqfile using recon, acqp, method and visu_pars files
% Ben Duffy Stanford University 2015

% files = dir;
% look at files.name to indices to folders
% list = {files(indstart:indend).name}'
% mkc_bd_readbruker(list)

% Added read header functions from aedes program written by Juha-Pekka Niskanen <Juha-Pekka.Niskanen@uef.fi>
% mkc 8/1/2018


inputfolderlist = cellstr(inputfolderlist);

cellfun(@readfolder,inputfolderlist)

end

function readfolder(foldername)

folderpath = [foldername];
if isunix; slash = '/'; elseif ispc; slash = '\'; end
pdatapath = [foldername slash 'pdata' slash '1' slash];
current = pwd;
if ispc
    processpath = [current(1:3) 'Data_Processed' slash 'fMRI'];
elseif isunix
    processpath = [current(1:end-4) 'Data_Processed' slash 'fMRI'];
end
fileoutputname = foldername
n_discard = 0; 


% Read Header files -------------------------

% Read ACQP file
acqp_file = [folderpath slash 'acqp'];
if exist(acqp_file,'file')~=2
    msg = sprintf('Cannot find ACQP file in "%s".',folderpath);
    return
end
try
    hdr.acqp = aedes_readjcamp(acqp_file);
catch
    hdr = [];
    msg = sprintf('Error while reading "%s". The error was "%s".',...
        acqp_file,lasterr);
    return
end

% Read Method file
method_file = [folderpath slash 'method'];
if exist(method_file,'file')~=2
    msg = sprintf('Cannot find METHOD file in "%s".',folderpath);
    return
end
try
    hdr.method = aedes_readjcamp(method_file);
catch
    hdr = [];
    msg = sprintf('Error while reading "%s". The error was "%s".',...
        method_file,lasterr);
    return
end

% Read Reco file
reco_file = [pdatapath slash 'reco'];
if exist(reco_file,'file')~=2
    msg = sprintf('Cannot find RECO file in "%s".',pdatapath);
    return
end
try
    hdr.reco = aedes_readjcamp(reco_file);
catch
    hdr = [];
    msg = sprintf('Error while reading "%s". The error was "%s".',...
        reco_file,lasterr);
    return
end

% Read visu_pars file
visu_pars_file = [pdatapath slash 'visu_pars'];
if exist(visu_pars_file,'file')~=2
    msg = sprintf('Cannot find VISU_PARS file in "%s".',pdatapath);
    return
end
try
    hdr.visu_pars = aedes_readjcamp(visu_pars_file);
catch
    hdr = [];
    msg = sprintf('Error while reading "%s". The error was "%s".',...
        visu_pars,lasterr);
    return
end

% Skip set up scans
if strcmp(hdr.acqp.PULPROG, 'FLASH.ppg') || strcmp(hdr.acqp.ACQ_protocol_name , '2_Localized_shim') || strcmp(hdr.acqp.PULPROG, 'STEAM.ppg')
    return
end

% Get parameters. Can actually get rid of all of this...
reco_ft1 = hdr.reco.RECO_ft_size(1);
reco_ft2 = hdr.reco.RECO_ft_size(2);
nsl = hdr.reco.RecoObjectsPerRepetition;
fov1 = hdr.reco.RECO_fov(1);
fov2 = hdr.reco.RECO_fov(2);
slthk = hdr.acqp.ACQ_slice_thick;

% Scale up so that SPM doesn't go crazy
vox_size1 = (fov1./reco_ft1)*100;
vox_size2 = (fov2./reco_ft2)*100;
vox_size3 = slthk*10;


fid = fopen([pdatapath '2dseq']);
if (fid == -1)
    return;
end
data = fread(fid,'uint16');
fclose(fid);

%Check to see what the appropriate params are for the conversion to .nii

if strcmp(hdr.acqp.ACQ_protocol_name , 'F_GE_EPI_1Rep_adj') ||...
        contains(hdr.acqp.ACQ_protocol_name,'EPI')
    
    % find number of reps and if short test scan then skip
    nreps = hdr.reco.RecoNumRepetitions;
    if nreps < 60
        return;
    end
    
    breaks = strfind(foldername,slash);
    if ispc
        datefolder = [processpath slash foldername(breaks(3)+1:breaks(3)+8) slash];
        resultfolder = [datefolder foldername(breaks(3)+17:breaks(4)-1) slash foldername(breaks(end):end) 'f'];
    elseif isunix
        datefolder = [processpath slash foldername(breaks(6)+1:breaks(6)+8) slash];
        resultfolder = [datefolder foldername(breaks(6)+17:breaks(7)-1) slash foldername(breaks(end):end) 'f'];
    end
    if ~exist(resultfolder,'dir'); mkdir(resultfolder); end
    
    fileoutputname = [resultfolder foldername(breaks(end):end) '_EPI'];
    
    display(sprintf('converting %s to .nii',foldername))
    data1 = reshape(data,reco_ft1,reco_ft2,nsl,nreps);
    for i = 1:size(data1,4)
        for j = 1:size(data1,3)
            data1(:,:,j,i) = fliplr(data1(:,:,j,i));
        end
    end
    
    if ndims(data1)==4
        data2 = data1(:,:,:,n_discard+1:end); %discard first n scans!!!!
        display(sprintf('warning: first %d (dummy) scans discarded',n_discard))
        nii = make_nii(data2,[vox_size1 vox_size2 vox_size3]);
    else
        nii = make_nii(data1,[vox_size1 vox_size2 vox_size3]);
    end
    save_nii(nii,[fileoutputname '.nii']);
    % sdi4(data2);
    
    
elseif strcmp(hdr.acqp.PULPROG , 'RARE.ppg')
    
   breaks = strfind(foldername,slash);
    if ispc
        datefolder = [processpath slash foldername(breaks(3)+1:breaks(3)+8) slash];
        resultfolder = [datefolder foldername(breaks(3)+17:breaks(4)-1) slash foldername(breaks(end):end) 'a'];
    elseif isunix
        datefolder = [processpath slash foldername(breaks(6)+1:breaks(6)+8) slash];
        resultfolder = [datefolder foldername(breaks(6)+17:breaks(7)-1) slash foldername(breaks(end):end) 'a'];
    end
    if ~exist(resultfolder,'dir'); mkdir(resultfolder); end
    
    data1 = reshape(data,reco_ft1,reco_ft2,nsl);
    for j = 1:size(data1,3)
        data1(:,:,j) = fliplr(data1(:,:,j));
    end
    nii = make_nii(data1,[vox_size1 vox_size2 vox_size3]); 
    fileoutputname = [resultfolder foldername(breaks(end):end) '_anatomy'];
    if reco_ft1 > 128
        fileoutputname = [resultfolder foldername(breaks(end):end) '_anatomy'];
    else
        fileoutputname = [resultfolder foldername(breaks(end):end) '_lowres_anatomy'];
    end  
    save_nii(nii,[fileoutputname '.nii']);
    
    
elseif strcmp(hdr.acqp.ACQ_protocol_name , 'fs_T2StarMap')
    
    % Have to load the map reconned in Paravision
    pdatapath2 = [pwd '\' foldername '\pdata\2\'];
    fid = fopen([pdatapath2 '2dseq']);
    if (fid == -1)
        return;
    end
    data = fread(fid,'int32','l');
    fclose(fid);
    
    % Have to load the corresponding params
    visu_pars_file = [pdatapath2,'\visu_pars'];
    hdr.visu_pars = aedes_readjcamp(visu_pars_file);
    data1 = reshape(data,reco_ft1,reco_ft2,hdr.visu_pars.VisuCoreFrameCount);
    
    for k=1:hdr.visu_pars.VisuCoreFrameCount
        data2(:,:,k)=data1(:,:,k)*hdr.visu_pars.VisuCoreDataSlope(k);
    end
    data3=reshape(data2,size(data2,1),size(data2,2),5,hdr.visu_pars.VisuCoreFrameCount/5);
    fileoutputname = [foldername '_T2s'];
    nii = make_nii(data3,[vox_size1 vox_size2 vox_size3]);
    save_nii(nii,[fileoutputname '.nii']);
    
elseif strcmp(hdr.acqp.ACQ_protocol_name , 'fs_B0Map')
   
    % Note that b0 map is int32 and not int16 like the other scans
    fid = fopen([pdatapath '2dseq']);
    if (fid == -1)
        return;
    end
    data = fread(fid,'int32','l');
    fclose(fid);
    data1 = (reshape(data,reco_ft1,reco_ft2,hdr.acqp.ACQ_size(3)))*hdr.visu_pars.VisuCoreDataSlope;
    fileoutputname = [foldername '_B0'];
    
    % For analyze format, which is necessary to use the FieldMap toolbox in
    % SPM
    
    %     nii = make_ana(data1,[vox_size1 vox_size2 vox_size3]);
    %     save_untouch_nii(nii,[fileoutputname '.nii']);
    % for nii format
    nii = make_nii(data1,[vox_size1 vox_size2 vox_size3]);
    save_nii(nii,[fileoutputname '.nii']);
    
else
    return;
end

save([fileoutputname '_params.mat'], 'hdr');

end

%% Reads in parameters from recon, acqp, method and visu_pars files and sticks them in a hdr file
% Juha-Pekka Niskanen <Juha-Pekka.Niskanen@uef.fi>

function jdx = aedes_readjcamp(filename)
% AEDES_READJCAMP - Read JCAMP DX format files (Bruker parameter files)
%
%
% Synopsis:
%       jdx=aedes_readjcamp(filename)
%
% Description:
%       The function reads the JCAMP DX files and returns a
%       structure with parameters as structure fields. The input
%       argument is a string containing the full path to the file.
%
% Examples:
%       jdx=aedes_readjcamp('C:\path\to\jcamp_dx_file')
%
% See also:
%       AEDES_READBRUKER, AEDES_DATA_READ, AEDES

jdx = [];

% Prompt for a file if not given as an input argument
if nargin == 0
    [fn,fp] = uigetfile({'*.*','All Files (*.*)'},'Open a JCAMP DX file');
    if isequal(fn,0)
        return
    end
    filename = [fp,fn];
elseif nargin > 1
    error('Too many input arguments.');
end

% Open the file for reading
fid = fopen(filename,'r');
if fid < 0
    error('Could not open file "%s" for reading.',filename);
end

% Check that the file is a JCAMP DX file
str = fread(fid,20,'char=>char');
if isempty(regexp(str.','^\s*##TITLE'))
    fclose(fid);
    error('File "%s" is not a valid JCAMP DX file.',filename)
end
fseek(fid,0,-1); % Rewind file

C = fread(fid,inf,'char');
fclose(fid);

% Remove carriage returns
C(C==13)=[];

% Convert to string
C = char(C.');

% Remove comment lines
C = regexprep(C,'\$\$([^\n]*)\n','');

% Remove unnecessary line breaks
f = @l_RemoveLineBreaks;
C=regexprep(C,'^(\s*[^#].*?)(?=\n\s*#)','${f($1)}','lineanchors');
C=regexprep(C,'(\([^\)]+?)\n(.*?\))','${f([$1,$2])}','lineanchors');
CC = regexp(C,'\s*##','split');
CC(1)=[];

% Parse the file line-by-line
for ii=1:length(CC)
    
    str = CC{ii};
    if strncmp(str,'END=',4)
        continue
    end
    
    % The commented regexp sometimes fails with long strings...
    %param = regexp(str,'^(.*)=','tokens','once');
    ind = find(str==61); % Find '=' chars...
    if isempty(ind)
        param='';
    else
        param=str(1:ind(1)-1);
    end
    %param = strrep(param{1},'$','');
    param = strrep(param,'$','');
    param = l_CheckParameter(param);
    
    if any(str==sprintf('\n'))
        % Get size
        sz = regexp(str,'=\s*\((.*?)\)\s*\n','tokens','once');
        sz = str2num(['[',sz{1},']']);
        
        % Parse value
        value = regexp(str,'\n(.*)$','tokens','once');
        value = value{1};
        value = l_CheckValue(value,sz);
    else
        value = regexp(str,'=\s*(.*)','tokens','once');
        value = value{1};
        value = l_CheckValue(value);
    end
    
    % Add to structure
    jdx.(param) = value;
    
end
end






% ==========================
% - Subfunctions -
% ==========================

% - Remove linebreaks
function out = l_RemoveLineBreaks(str)

out = strrep(str,sprintf('\n'),'');

end

% - Check parameter value --------------------------
function out = l_CheckValue(val,sz)

if nargin == 1
    sz = 0;
end

% Remove insignificant whitespace
val = strtrim(val);

if isempty(val)
    out = val;
    return
end

% Handle strings and string lists
if val(1) == '<' && val(end) == '>'
    val(val=='<')='''';
    val(val=='>')='''';
    out = eval(['{',val,'}']);
    if length(out) == 1
        out = out{1};
    end
    return
end

% Handle cell matrices
if val(1) == '(' && val(end) == ')'
    nRows = length(find(val==')'));
    
    % Nested tables are not supported. This is a workaround for nested tables
    % and everything is read in a single lined table...
    if nRows ~= sz && sz>0
        nRows=sz;
    end
    
    val(1) = '';
    val(end) = '';
    val(val=='(')='';
    val(val==')')=',';
    val(val=='<')='';
    val(val=='>')='';
    
    % Split using the commas
    val_split = regexp(val,',\s+','split');
    val_out = cell(size(val_split));
    
    % Try to convert to numbers
    for ii = 1:length(val_split)
        num = str2double(val_split{ii});
        if isnan(num)
            val_out{ii} = val_split{ii};
        else
            val_out{ii} = num;
        end
    end
    
    
    out = reshape(val_out,[],nRows).';
    return
end

% Check if the string contains only numbers before tryin to convert to a
% number. str2num uses eval command and if the string matches to a
% function name strange things can happen...
tmp2 = regexp(val,'[^\d\.\seE-+]');
if ~isempty(tmp2)
    out = val;
    return
end

% Convert value to numeric if possible
tmp = str2num(val);
if ~isempty(tmp) && isreal(tmp)
    if length(sz)>1
        tmp = reshape(tmp,sz(2),sz(1),[]);
        tmp = permute(tmp,[2 1 3]);
    end
    out = tmp;
    return
end

out = val;
end

% - Check parameter strings -------------------------
function out = l_CheckParameter(param)

alphabets = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ';
numbers = '1234567890';

% Remove insignificant whitespace
param = strtrim(param);

if isempty(param)
    out = 'EMPTY_PARAM';
    return
end

% Check parameter starts with a valid structure field character
if ~any(param(1)==alphabets)
    param = ['PAR_',param];
end

% Check that the parameter string does not contain any illegal characters
% (for Matlab structure fields)
ind = ~ismember(param,[alphabets,numbers,'_']);
if any(ind)
    param(ind) = '_';
end

out = param;
end
