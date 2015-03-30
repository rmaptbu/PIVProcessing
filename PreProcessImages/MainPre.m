clear;clc;
addpath('C:\Users\jsherwoo\Desktop\MatlabProcessing\PreProcessImages')          % path of sub-functions etc.
basepathS='C:\Users\jsherwoo\Desktop\MatlabProcessing\SampleCase\Raw\Str\';     % path for raw strobe images
basepathL='C:\Users\jsherwoo\Desktop\MatlabProcessing\SampleCase\Raw\Las\';     % path for laser images
writepath='C:\Users\jsherwoo\Desktop\MatlabProcessing\SampleCase\Processed';   % path to save processed data
cd(basepathS)
cases=struct2cell(dir);                             % get a list of the file names
cases=cases(1,3:end);
for i=1:length(cases)                               % for each case
    disp(['Case ', num2str(i)])
    basepathScase=[basepathS,char(cases(i)),'\'];   % Build list of strobe images
    cd(basepathScase)
    filesS=struct2cell(dir);
    filesS=filesS(1,3:end);
    basepathLcase=[basepathL,char(cases(i)),'\'];   % Build list of laser images
    cd(basepathLcase) 
    filesL=struct2cell(dir);
    filesL=filesL(1,3:end);
    disp('Locate and Mask')
    [Columns,Rows,ImMeanF,Ga,Angle,Change,LasImMin]=Preprocess(basepathScase,filesS,basepathLcase,filesL);
    disp('Write files')
    % save mean image
    if isdir([writepath,'\Mean\'])==0;mkdir([writepath,'\Mean\']);end
    ImMeanFR=fliplr(imresize(ImMeanF,84/77)); % resize to diameter of 84 pixels
    imwrite(uint16(ImMeanFR),[writepath,'\Mean\',char(cases(i)),'_mean.tif'])
    
    % save individual images
    if isdir([writepath,'\Las\'])==0;mkdir([writepath,'\Las\']);end
    if isdir([writepath,'\Str\'])==0;mkdir([writepath,'\Str\']);end

    WritePathL=[writepath,'\Las\',char(cases(i)),'\'];mkdir(WritePathL);
    WritePathS=[writepath,'\Str\',char(cases(i)),'\'];mkdir(WritePathS);
    % write files
    for j=1:length(filesS)
        disp(['Write image ', num2str(j)])
        ImS=double(imread([basepathScase,char(filesS(j))])-32768);  % import images
        ImL=double(imread([basepathLcase,char(filesL(j))])-32768);
        H=size(ImS,1)/2; 
        ImSA=ImS(1:H,:);ImSB=ImS(H+1:end,:);                    % split images in two
        ImLA=ImL(1:H,:);ImLB=ImL(H+1:end,:); 
        %% correct images
        ImSA=ImSA./Ga;ImSB=ImSB./Ga;                % correct strobe for illumination
        ImLA=ImLA-LasImMin;ImLB=ImLB-LasImMin;      % subtract minimium imagefrom laser
        % rotate images, crop excess 
        ImSA=imrotate(ImSA,Angle,'bicubic');ImSA=ImSA(Change+1:end-Change,Change+1:end-Change);
        ImSB=imrotate(ImSB,Angle,'bicubic');ImSB=ImSB(Change+1:end-Change,Change+1:end-Change);
        ImLA=imrotate(ImLA,Angle,'bicubic');ImLA=ImLA(Change+1:end-Change,Change+1:end-Change);
        ImLB=imrotate(ImLB,Angle,'bicubic');ImLB=ImLB(Change+1:end-Change,Change+1:end-Change);
        %% Crop to region
        ImSA=ImSA(Columns,Rows);ImSB=ImSB(Columns,Rows);
        ImLA=ImLA(Columns,Rows);ImLB=ImLB(Columns,Rows);
        %% resize
        ImSAr=imresize(ImSA,84/77);ImSBr=imresize(ImSB,84/77);
        ImLAr=imresize(ImLA,84/77);ImLBr=imresize(ImLB,84/77);
        
        ImSwrite=fliplr([ImSAr;ImSBr]);
        ImLwrite=fliplr([ImLAr;ImLBr]);
        
        imwrite(uint16(ImLwrite),[WritePathL,char(filesL(j))]);
        imwrite(uint16(ImSwrite),[WritePathS,char(filesS(j))]);
    end
end

        
    
 