clear;clc;
% path of sub-functions etc.
addpath('/Users/Thore/Documents/MATLAB/PIVProcessing/PreProcessImages/');
% path for raw strobe images
basepath='/Users/Thore/Documents/PIVLabData/oct10_10muSpheres_400muTube/Str/';     
% path to save processed data
writepath='/Users/Thore/Documents/PIVLabData/oct10_10muSpheres_400muTube/Processed/';   
cd(basepath)
cases=struct2cell(dir);                             % get a list of the file names
cases=cases(1,3:end);
for i=1:length(cases)                               % for each case
    disp(['Case ', num2str(i)])
    % Build list of strobe images
    basepathScase=[basepath,char(cases(i)),'/'];   
    cd(basepathScase)
    filesS=struct2cell(dir);
    filesS=filesS(1,3:end);
    disp('Locate and Mask')
    [Columns,Rows,ImMeanF,Ga,Angle,Change,ImMax,IntCorrection]=Preprocess(basepathScase,filesS);
    disp('Write files')
    % save mean image
    if isdir([writepath,'/Mean/'])==0;mkdir([writepath,'/Mean/']);end
    ImMeanFR=fliplr(imresize(ImMeanF,400/401));
    ImMaxFR=fliplr(imresize(ImMax,400/401));
    imwrite(uint16(ImMeanFR),[writepath,'/Mean/',char(cases(i)),'_mean.tif'])
    imwrite(uint16(ImMax),[writepath,'/Mean/',char(cases(i)),'_max.tif'])
    
    % save individual images
    if isdir([writepath,'/Str/'])==0;mkdir([writepath,'/Str/']);end

    WritePathS=[writepath,'/Str/',char(cases(i)),'/'];mkdir(WritePathS);
    % write files
    for j=1:length(filesS)
        disp(['Write image ', num2str(j)])
        ImS=double(imread([basepathScase,char(filesS(j))])-32768);  % import images
        %% correct images
        ImS=ImS./Ga;         % correct strobe for illumination
        if mod(j,2)==0; ImS=ImS*IntCorrection;end            
        % rotate images, crop excess 
        ImS=imrotate(ImS,Angle,'bicubic','crop');
        ImS=ImS(Change+1:end-Change-1,Change+1:end-Change);
        % Crop to region
        ImS=ImS(Columns,Rows);
        % resize
        ImSr=imresize(ImS,400/401);        
        ImSwrite=fliplr([ImSr]);      
        imwrite(uint16(ImSwrite),[WritePathS,char(filesS(j))]);
    end
end

% for i=1:length(cases)                               % for each case
%     disp(['Case ', num2str(i)])
%     basepathScase=[basepathS,char(cases(i)),'\cropped\'];   % Build list of strobe images
%     cd(basepathScase)
%     filesS=struct2cell(dir);
%     filesS=filesS(1,3:end);
% for j=1:length(filesS);ImS=double(imread([basepathScase,char(filesS(j))]));
%ImSr=imresize(ImS,800/size(ImS,1));
%imwrite(uint16(ImSr),[basepathScase,'c',char(filesS(j))]);
% end;
% end
% 
%         
    
 