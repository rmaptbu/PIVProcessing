clear;clc;
addpath('G:\New\Image Preprocessing\PreProcessImages')
basepathS='G:\New\Images\Calibration\A\Raw\Str\';
writepath='G:\New\Images\Calibration\A\PreProcess\';
cd(basepathS)
cases=struct2cell(dir);
cases=cases(1,3:end);
for i=1:length(cases)   % for each case
    disp(['Case ', num2str(i)])
    basepathScase=[basepathS,char(cases(i)),'\'];
    cd(basepathScase)
    filesS=struct2cell(dir);
    filesS=filesS(1,3:end);
    disp('Locate and Mask')
    [Columns,Rows,ImMeanF,Ga,Angle,Change,ImMadF]=PreprocessCalib(basepathScase,filesS);
    disp('Write files')
    % save mean image
    if isdir([writepath,'Mean\'])==0;mkdir([writepath,'Mean\']);end
    ImMeanFR=fliplr(imresize(ImMeanF,84/77)); % resize to diameter of 84 pixels
    imwrite(uint16(ImMeanFR),[writepath,'\Mean\',char(cases(i)),'_mean.tif'])
    if ~isdir([writepath,'Mad\']);mkdir([writepath,'Mad\']);end
    ImMadFR=fliplr(imresize(ImMadF,84/77)); % resize to diameter of 84 pixels
    imwrite(uint16(ImMadFR),[writepath,'\Mad\',char(cases(i)),'_mad.tif'])

    % save individual images
    if isdir([writepath,'\Str\'])==0;mkdir([writepath,'\Str\']);end
    WritePathS=[writepath,'\Str\',char(cases(i)),'\'];mkdir(WritePathS);
    % write files
    for j=1:length(filesS)
        disp(['Write image ', num2str(j)])
        ImS=double(imread([basepathScase,char(filesS(j))])-32768);  % import images
        H=size(ImS,1)/2; 
        ImSA=ImS(1:H,:);ImSB=ImS(H+1:end,:);                    % split images in two
        %% correct images
        ImSA=ImSA./Ga;ImSB=ImSB./Ga;                % correct strobe for illumination
        % rotate images, crop excess 
        ImSA=imrotate(ImSA,Angle,'bicubic');ImSA=ImSA(Change+1:end-Change,Change+1:end-Change);
        ImSB=imrotate(ImSB,Angle,'bicubic');ImSB=ImSB(Change+1:end-Change,Change+1:end-Change);
        %% Crop to region
        ImSA=ImSA(Columns,Rows);ImSB=ImSB(Columns,Rows);
        %% resize
        ImSAr=imresize(ImSA,84/77);ImSBr=imresize(ImSB,84/77);
        ImSwrite=fliplr([ImSAr;ImSBr]);
        % save
        imwrite(uint16(ImSwrite),[WritePathS,char(filesS(j))]);
    end
end

        
    
 