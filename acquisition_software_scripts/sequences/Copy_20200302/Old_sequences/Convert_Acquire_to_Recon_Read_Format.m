% save the RFData collected at high frame rate

folder = 'Data/';
files = dir(folder);
filename = {files.name};

for i=3:size(filename,2)
    if ~isempty(strfind(filename{i},'P'))
        pFile=[folder filename{i}];
        dataFile=[folder filename{i}(1:end-6) '.mat'];
    
        datadir=dataFile(1:end-4);
        if ~exist(datadir,'dir')
            mkdir(datadir);
        end
            
        RFtmp = load(dataFile);
        load(pFile);

        movefile(dataFile,[datadir dataFile(5:end)]);
        movefile(pFile,[datadir pFile(5:end)]);

        disp ('Reconfiguring RF Data -- please wait!'), disp(' ')
        names = fieldnames(RFtmp);
        RcvData=cell(1,1);

        NumIter=ceil(size(names,1)/50);
        frame=1;
        for j=1:NumIter

            for k = 1:50
                if frame<=size(names,1)
                    RcvData{1,1}(:,:,k)=RFtmp.(names{frame});
                    frame = frame + 1;
                end
            end

            save ([datadir '/RFdataFlashAnglesHFR_' num2str(NumIter)], 'RcvData', 'P')
        end
    end
end
