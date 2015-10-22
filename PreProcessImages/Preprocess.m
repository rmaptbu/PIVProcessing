function[Columns,Rows,ImMeanFinal,Ga,Angle,Change,ImMax,IntCorrection]=Preprocess(basepathScase,filesS)
%% get mean and std images from strobe
Im=imread([basepathScase,char(filesS(1))]);     % read initial image for paramneters
H=size(Im,1);                                 % find height of single image
for l=1:length(filesS)                           
    Im=single(imread([basepathScase,char(filesS(l))])-32768); % correct to 12 bit
    ImStack(:,:,l)=Im;                      % get dif image
end
disp('Std image')
ImStd=std(ImStack,[],3);
disp('Mean image')
ImMean=mean(ImStack,3);
ImMax=max(ImStack,[],3);
ImMeanA= mean(mean(mean(ImStack(:,:,1:2:end-1))));
ImMeanB= mean(mean(mean(ImStack(:,:,2:2:end))));
IntCorrection=ImMeanA/ImMeanB;
clear ImStack

%% Fix Illumination for strobe images
element=[ 0 1 0; 1 1 1; 0 1 0];                     % mask out the channel
ImStdBW=im2bw(uint16(ImStd),8E-4);
ImStdBWD=imdilate(ImStdBW,element);for i=1:4;ImStdBWD=imdilate(ImStdBWD,element);end;
ImStdBWDE=imerode(ImStdBWD,element);for i=1:15;ImStdBWDE=imerode(ImStdBWDE,element);end;
disp('Fix Illumination')
Ga=FixStrobeIllumination(ImMean,ImStdBWDE); 
ImMeanF=ImMean./Ga;
disp('Fix Location')
%% Fix rotation
Angle=-39.55;
Change=tan(4.17/180*pi)*1024;
Change=ceil(abs(Change))+3;                 % this is the amount extra added to the image which should be trimmed off
                                            % trims an extra 3 pixels for safety                
ImMeanR=imrotate(ImMeanF,Angle,'bicubic','crop');ImMeanR=ImMeanR(Change+1:end-Change-1,Change+1:end-Change);

%% Crop Mean Image 
Rows=100:size(ImMeanR,2);Columns=280:680;
ImMeanFinal=ImMeanR(Columns,Rows);
imagesc(ImMeanFinal)

