%function[Columns,Rows,ImMeanFinal,Ga,Angle,Change,ImMin]=Preprocess(basepathScase,filesS,basepathLcase,filesL)
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
ImMin=min(ImStack,[],3);
clear ImStack
%% get min image for laser
% for l=1:length(filesL)                          
%     Im=single(imread([basepathLcase,char(filesL(l))])-32768); % correct to 12 bit
%     ImA=Im(1:H,:);ImB=Im(H+1:end,:);                    % Split Image       
%     ImStack(:,:,2*l-1)=ImA;                      % get dif image
%     ImStack(:,:,2*l)=ImB;  
% end
% disp('Laser Min image')
% LasImMin=min(ImStack,[],3);
% clear ImStack
%% Fix rotation
% ImStdBWE=imerode(ImStdBWD,element);
% for l=1:20;ImStdBWE=imerode(ImStdBWE,[1 1 1; 1 1 1; 1 1 1]);end
% ImStdLine=bwmorph(ImStdBWE,'skel',inf);     % Reduce to a line
% 
% [~,LHS]=max(ImStdLine(:,64));               % get line centre on left side
% [~,RHS]=max(ImStdLine(:,1281));             % get line centre on left side
% 
% Change=ceil(RHS-LHS);                       % find difference in long angle between sides
% Angle=atan2(Change,1216)*180/pi;            % calculate angle for rotation
Angle=4.71;
Change=tan(4.71/180*pi)*1024;
Change=ceil(abs(Change))+3;                 % this is the amount extra added to the image which should be trimmed off
                                            % trims an extra 3 pixels for safety                
ImMeanR=imrotate(ImMean,Angle,'bicubic','crop');ImMeanR=ImMeanR(Change+31:end-Change-31,Change+1:end-Change);
ImStdR=imrotate(ImStd,Angle,'bicubic','crop');ImStdR=ImStdR(Change+31:end-Change-31,Change+1:end-Change);
%LasImMinR=imrotate(LasImMin,Angle,'bicubic');LasImMinR=LasImMinR(Change+1:end-Change,Change+1:end-Change);
%% Fix Illumination for strobe images
element=[ 0 1 0; 1 1 1; 0 1 0];                     % mask out the channel
ImStdBW=im2bw(uint16(ImStdR),4E-4);
ImStdBWD=imdilate(ImStdBW,element);
for i=1:10;ImStdBWD=imdilate(ImStdBWD,element);end   %
disp('Fix Illumination')
Ga=FixStrobeIllumination(ImMeanR,ImStdBWD); 
ImMeanF=ImMeanR./Ga;
disp('Fix Location')
%% find parent centre

% f=fspecial('prewitt');
% edge=imfilter(ImMeanR,f,'replicate');
% prof=mean(edge(:,[1:300,600:800,1050:end]),2);
% [midP,wp]=FindMidline(prof);
% midP=floor(midP);
%% find branch centre by rotating filter
% rotedge=imfilter(ImMeanR,f','replicate');
% leftprof=-mean(rotedge(1:700,1:700),1);
% [midL,wl]=FindMidline(leftprof);
% rightprof=-mean(rotedge(1:700,701:end),1);
% [midR,wr]=FindMidline(rightprof);
% midR=midR+700;
% D=77; % probably slightly less, but this seems closest
% imagesc(ImMeanR);axis image;colormap(gray);
% line([0 1344],[midP midP]);line([0 1344],[midP+0.5*D midP+0.5*D]);line([0 1344],[midP-0.5*D midP-0.5*D]);
% line([midL midL],[0 1024]);line([midL+0.5*D midL+0.5*D],[0 1024]);line([midL-0.5*D midL-0.5*D],[0 1024]);
% line([midR midR],[0 1024]);line([midR+0.5*D midR+0.5*D],[0 1024]);line([midR-0.5*D midR-0.5*D],[0 1024]);

% centre of image is 
% CentreX=floor((midL+midR)/2); 
% ParentY=midP;

%% get regions 
% xl=CentreX-8*D+1;     % gives 5 diameters from channel centre
% xr=CentreX+8*D;     % gives 5 diameters from channel centre
% yu=ParentY-10*D+1;     % gives daughter branches of 10 diameters
% yl=ParentY+D;       % gives half a diameter extra at bottom of image
Rows=1:size(ImMeanR,1);Columns=1:size(ImMeanR,2);
%% Crop Mean Image 
ImMeanFinal=ImMeanR(Columns,Rows);

