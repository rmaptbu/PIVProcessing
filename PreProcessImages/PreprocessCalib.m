function[Columns,Rows,ImMeanFinal,Ga,Angle,Change,ImMadFinal]=PreprocessCalib(basepathScase,filesS)
%% get mean and std images from strobe
Im=imread([basepathScase,char(filesS(1))]);     % read initial image for paramneters
H=size(Im,1)/2;                                 % find height of single image
for l=1:length(filesS)                           
    Im=single(imread([basepathScase,char(filesS(l))])-32768); % correct to 12 bit
    ImA=Im(1:H,:);ImB=Im(H+1:end,:);                    % Split Image
    ImStack(:,:,2*l-1)=ImA;                      % get dif image
    ImStack(:,:,2*l)=ImB;  
end
disp('Std image')
ImStd=std(ImStack,[],3);
ImMad=mad(ImStack,[],3);
disp('Mean image')
ImMean=mean(ImStack,3);
clear ImStack
%% Fix Illumination for strobe images
element=[ 0 1 0; 1 1 1; 0 1 0];                     % mask out the channel
ImStdBW=im2bw(uint16(ImStd),100/(2^16));
ImStdBWD=imdilate(ImStdBW,element);
for i=1:4;ImStdBWD=imdilate(ImStdBWD,element);end   %
disp('Fix Illumination')
Ga=FixStrobeIllumination(ImMean,ImStdBWD); 
ImMeanF=ImMean./Ga;
ImMadF=ImMad./Ga;
disp('Fix Location')
%% Fix rotation
ImStdBWE=imerode(ImStdBWD,element);
for l=1:10;ImStdBWE=imerode(ImStdBWE,[1 1 1; 1 1 1; 1 1 1]);end
ImStdLine=bwmorph(ImStdBWE,'skel',inf);     % Reduce to a line

[~,LHS]=max(ImStdLine(:,64));               % get line centre on left side
[~,RHS]=max(ImStdLine(:,1281));             % get line centre on left side

Change=ceil(RHS-LHS);                       % find difference in long angle between sides
Angle=atan2(Change,1216)*180/pi;            % calculate angle for rotation
Change=ceil(abs(Change))+3;                 % this is the amount extra added to the image which should be trimmed off
                                            % trims an extra 3 pixels for safety                
ImMeanR=imrotate(ImMeanF,Angle,'bicubic');ImMeanR=ImMeanR(Change+1:end-Change,Change+1:end-Change);
ImMadR=imrotate(ImMadF,Angle,'bicubic');ImMadR=ImMadR(Change+1:end-Change,Change+1:end-Change);
%% find parent centre
f=fspecial('prewitt');
edge=imfilter(ImMeanR,f,'replicate');
prof=mean(edge(:,[1:300,600:800,1050:end]),2);
[midP]=FindMidlineCalib(prof);
midP=floor(midP);
D=77; % probably slightly less, but this seems closest
imagesc(ImMeanR);axis image;colormap(gray);
line([0 1344],[midP midP]);line([0 1344],[midP+0.5*D midP+0.5*D]);line([0 1344],[midP-0.5*D midP-0.5*D]);

% centre of image is 
CentreX=ceil(size(ImMeanR,2)/2);
ParentY=midP;

%% get regions 
xl=CentreX-5*D+1;     % gives 5 diameters from channel centre
xr=CentreX+5*D;     % gives 5 diameters from channel centre
yu=ParentY-D+1;     % gives daughter branches of 10 diameters
yl=ParentY+D;       % gives half a diameter extra at bottom of image
Rows=xl:xr;Columns=yu:yl;
%% Crop Mean Image 
ImMeanFinal=ImMeanR(Columns,Rows);
ImMadFinal=ImMadR(Columns,Rows);

