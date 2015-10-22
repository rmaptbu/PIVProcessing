%% load files
clear;

addpath('/Users/Thore/Documents/MATLAB/PIVProcessing/');
basepath='/Users/Thore/Documents/PIVLabData/oct10_10muSpheres_400muTube/Processed/Results/';  
cd(basepath)
cases=struct2cell(dir);              % get a list of the file names
cases=cases(1,3:end);
names=[];

maps=[];
for i=1:length(cases)
    if strfind(cases{i},'.jvc')
        maps=cat(4,maps,load(cases{i}));
        name=cases{i}
        names=cat(1,names,name);
    end
end

maps=maps(:,3:4,:,:);
NvecX=133;%137
for i=1:NvecX:length(maps(:,1,1,1))
    maps(1:NvecX,:,(i+NvecX-1)/NvecX,:)=maps(i:i+NvecX-1,:,1,:);
end
maps=maps(1:NvecX,:,:,:);

%% finde mean speed -> maps(X,v,Y,file)
meanspeed=(maps(:,1,:,:).^2+maps(:,2,:,:).^2).^0.5;
meanspeed=squeeze(meanspeed)*0.63; %convert from pixel to microns
% figure;subplot(2,2,1:2)
%surf(meanspeed(:,:,1))
meanspeed=squeeze(mean(meanspeed,1));
dt=[4 8 12 16 20 40 50 60 80];
for i=1:length(dt)
meanspeed2(:,i)=meanspeed(:,i)./dt(i);
end


figure;plot(meanspeed);title('Mean displacement')
figure;plot(meanspeed2);title('Mean velocity')

 for i=1:length(dt); subplot(4,4,i);     
     plot(meanspeed(:,i), 'Color', [0, 0, 0]);
     S=names(i,:);
     title(S);
 end