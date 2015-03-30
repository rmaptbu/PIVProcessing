% Get Average Illumination
% for calib - 
        % get channel values and normalise by middle value, then multiply
        % by 4095
clear;clc;
D=84;
ParentX=1:4.5*D;            % parent branch
ParentY=9.5*D+1:10.5*D;
MidX=6.5*D+1:9.5*D;         % section in channel centre for quanitfying overall brightness of image
MidY=5.5*D+1:8.5*D;
Trim=6;
basepath=['G:\New\Images\A\PreProcess_Double\Str\';'G:\New\Images\N\PreProcess_Double\Str\'];
calibpath=['G:\New\Images\Calibration\A\PreProcess\Str\';'G:\New\Images\Calibration\N\PreProcess\Str\'];
PC=0.1;
% cc=0.95; % correction value - difference between outer channel and inner channel with no haematocrit.
% % Use histogram normalisation for translating data to be same as
% % calibration

for Agg=1:2
    %%  calibration images
    cd(calibpath(Agg,:));
    calib=struct2cell(dir);calib=calib(1,3:end);
    k=1;
    for i=1:length(calib);                                  %  get calibration images for H=0.25 cases
        if isempty(findstr(char(calib(i)),'25'))==0;
            Calib25(k)=calib(i);k=k+1;
        end
    end
    for i=1:numel(Calib25)                                  % for each of the three H=0.25 cases
        cd(char(Calib25(i)))
        images=struct2cell(dir);images=images(1,3:end);
        for j=1:length(images)
            I=single(imread(char(images(j))));                                          % import the image
            I=I(:,3*D+1:7*D);                                                           % select inner 4 diameters
            CalibOut=I([1:ceil(0.4*D),ceil(1.6*D):ceil(2.4*D),ceil(3.6*D):end],:);      % get the outer region of the image
            calibMid(j)=mean(CalibOut(:));                                              % outer brightness for calib images
            Itrim=I([D/2+Trim+1:3*D/2-Trim,5*D/2+Trim+1:7*D/2-Trim],:);                 % select the channel
            IHist(:,j,i)=Itrim(:); 
%             bar(hist(IHist(:,j,i),0:2:4095)/numel(IHist(:,j,i)));
        end
        cd ..
    end
 
    CalibMid=mean(calibMid(:));   
    IHistn=IHist/CalibMid*4095;         % Normalise calibration histogram so that channel outer region is 4095
     
    lowout(Agg)=prctile(IHistn(:),PC)/4095;       % use PC percentile for scaling
    highout(Agg)=coerce1(prctile(IHistn(:),100-PC)/4095,0,1);
%     highout=[1 1];
                      % average across the three calibration images
    %% now for individual cases - import individual images get histogram in parent branch
    cd(basepath(Agg,:));    
    Cases=struct2cell(dir);Cases=Cases(1,3:end);
    for i=1:length(Cases)
        cd(char(Cases(i)))
        images=struct2cell(dir);images=images(1,3:end);
        for j=1:length(images)
             [Agg Cases(i) j 'a']
            Im=imread(char(images(j)));   
            ImA=Im(1:size(Im,1)/2,:);
            ImB=Im(size(Im,1)/2+1:end,:);
            ImMidA=ImA(MidY,MidX);    
            ImMidB=ImB(MidY,MidX);  
            imMidA(j,i)=mean([ ImMidA(:); ImMidB(:)]);
            ParA=single(ImA(ParentY(Trim+1:end-Trim),ParentX));     % get parent branch
            ParB=single(ImB(ParentY(Trim+1:end-Trim),ParentX));   
            ParHist(:,j)=[ParA(:);ParB(:)];
        end
        
        ImMid=mean(imMidA(:,i));  
        ParHistn=ParHist/(ImMid)*4095;  % Normalise parent histogram so that channel outer region is 4095
%         subplot(312)
%         bar(0:2:4095,hist(ParHistn(:),0:2:4095)/numel(ParHistn(:)));
%              
        lowin=prctile(ParHistn(:),PC)/4095;                      % get lower 1%
        highin=coerce1(prctile(ParHistn(:),100-PC)/4095,0,1); 
        
        path=[basepath(Agg,1:end-4),'meanCor\',char(Cases(i)),'_mean_cor.tif'];
         for j=1:length(images)
             [Agg Cases(i) j 'b']
            Im=coerce2(single(imread(char(images(j))))/single(imMidA(j,i))*4095,0,4095); % normalise image
            ImA=Im(1:size(Im,1)/2,:);       % split image
            ImB=Im(size(Im,1)/2+1:end,:);
            ImAa=imadjust(ImA/4095,[lowin; highin],[lowout(Agg); highout(Agg)]);    % adjust both images
            ImBa=imadjust(ImB/4095,[lowin; highin],[lowout(Agg); highout(Agg)]);
            ImCA(:,:,j)=ImAa;%/max(ImAa(:))*4095;%./(max(ImAa(:)));   % re-normalising- dont do this
            ImCB(:,:,j)=ImBa;%/max(ImAa(:))*4095;%./(max(ImBa(:)));
         end
         %
%             a=Itrim(1:72,:)/CalibMid*4095;
%             b=ImA(ParentY(Trim+1:end-Trim),ParentX);
%             c=ImAa(ParentY(Trim+1:end-Trim),ParentX)*4095;    
% 
%         figure(1)
%           subplot(311);imagesc(a);caxis([0 4095])
%           subplot(312);imagesc(b);caxis([0 4095])
%           subplot(313);imagesc(c);caxis([0 4095])
% 
%         figure(2)
%           subplot(211); bar(0:2:4095,hist(a(:),0:2:4095)/numel(a));axis([0 4095 0 0.01])
%           subplot(212); bar(0:2:4095,hist(b(:),0:2:4095)/numel(b));axis([0 4095 0 0.01]);hold on
%           bar(0:2:4095,hist(c(:),0:2:4095)/numel(c),'r');axis([0 4095 0 0.01])

         ImC(:,:,1:size(ImCA,3))=ImCA;                  % recombine images
         ImC(:,:,size(ImCA,3)+1:2*size(ImCA,3))=ImCB;
         ImMean=uint16(mean(ImC,3)*4095);               % scale back up to 12-bit
         ImSem= uint16(std(ImC,[],3)./sqrt(size(ImC,3))*4095); 
%          ParA2=single(ImMean(ParentY(Trim+1:end-Trim),ParentX));     % get parent branch
%          ParHistC=ParA2(:);
%          
%          subplot(313)
%          bar(hist(ParHistC,0:2:4095)/numel(ParHistC));
%              
         
%          if Agg==2; ImMean=ImMean;end
         imwrite(ImMean,path)
         cd ..
    end
end

        
            