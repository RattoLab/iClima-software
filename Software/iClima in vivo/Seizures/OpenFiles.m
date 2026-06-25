function  OpenFiles(app)
% Open new data file according to the selected data type

% Check if any accessory window is selected and, if that is the case, closes it.

sel_format = get(app.dataFormat,'Value');
switch sel_format
    case 'National Instr'
        fmtSt = '*.dat';
    case 'EEG Laura'
        fmtSt = '*.eeg';
    case 'Meyer VEP'
        fmtSt = '*.asc';
    case 'SEP'
        fmtSt = '*.asc';
    case 'AxoP'
        fmtSt = '*.abf';
    case 'TRC import'
        fmtSt = '*.trc';
    case 'IMM MEA'
        fmtSt = '*.mat';
    case 'EDF format'
        fmtSt = '*.edf';
    case 'NewCastle'
        fmtSt = '*.mat';
    case 'Neuronexus'
        fmtSt = '.json';
    case 'Matlab'
        fmtSt = '.mat';
end
selectedFolder=uigetdir(app.dir_in);
folderContents=dir(selectedFolder);
 figure(app.ZebraMainFig)
for i = 1:length(folderContents)
    % Check if the file has the .abf extension
    [~,~,ext] = fileparts(folderContents(i).name);
    ext = ['*', ext];
    if strcmpi(ext, fmtSt)
       
      filePath = fullfile(selectedFolder, folderContents(i).name);
      fileInfo = {'name' 'date' 'bytes' 'isdir' 'datenum'};
        
        fileInfo = dir(filePath);
        timeNum = datestr(fileInfo.datenum);
        app.fileTime = timeNum (13:20);
        openFile(app,folderContents(i).name,[selectedFolder filesep])
    end
end
end

       