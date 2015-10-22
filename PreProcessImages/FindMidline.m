function[mid,w]=FindMidline(prof)
%% simplified code - this just find the two largest peaks and uses the centre
% there must have been a very repeatable pattern.
% find  positions of edges
% [sizes,edges]=findpeaks(double(abs(prof)),'MINPEAKHEIGHT',1000,'MINPEAKDISTANCE',50);   % find left edge
% lredges=edges(end-1:end);    % select two last peaks
% 
% ledge=min(lredges);     % select left edge
% redge=max(lredges);     % select right edge
% 
% mid=(ledge+redge)/2;
% % w=R-L;

%% applied code
% [sizes,edges]=findpeaks(double(abs(prof)),'MINPEAKHEIGHT',1000,'MINPEAKDISTANCE',50);   % find left edge
% % [m]=sort(sizes);      % put in order of size - this is wrong :(
% lredges=edges(end-1:end);    % select two largest peaks
% 
% ledge=min(lredges);     % select left edge
% redge=max(lredges);     % select right edge
% % ledge=edges; CalibN H10 high - peak finding did not work
% % redge=319 % 
% Buffer=20;                    % this bit is ignored....
% [~,L]=min(prof(ledge-Buffer:ledge+Buffer));L=L+ledge-1-Buffer;
% [~,R]=min(prof(redge-Buffer:redge+Buffer));R=R+redge-1-Buffer;
% 
% mid=(ledge+redge)/2;          % should have been L and R
% w=R-L;

%% corrected applied code
[sizes,edges]=findpeaks(double(abs(prof)),'MINPEAKHEIGHT',100,'MINPEAKDISTANCE',50);   % find left edge
[m,I]=sort(sizes);              % put in order of size 
lredges=edges(I(end-1:end));    % select two largest peaks

ledge=min(lredges);     % select left edge
redge=max(lredges);     % select right edge

Buffer=20;
[~,L]=min(prof(ledge-Buffer:ledge+Buffer));L=L+ledge-1-Buffer;
[~,R]=min(prof(redge-Buffer:redge+Buffer));R=R+redge-1-Buffer;

mid=(L+R)/2;
w=R-L;















% for parent branch
% subplot(1,2,1);imagesc(im);line([0 1344],[mid mid]);colormap(gray)
% line([0 1344],[L L ]);line([0 1344],[R R]);
% D=80;
% subplot(1,2,2);imagesc(ImMeanR);line([0 1344],[mid mid]);colormap(gray)
% line([0 1344],[mid+0.5*D mid+0.5*D]);line([0 1344],[mid-0.5*D mid-0.5*D]);

% for left branch
% subplot(132);imagesc(im);colormap(gray)
% line([mid mid],[0 1024]);line([L L],[0 1024]);line([R R],[0 1024]);
% D=80;
% subplot(133);imagesc(ImMeanR);
% line([mid mid],[0 1024]);line([mid+0.5*D mid+0.5*D],[0 1024]);line([mid-0.5*D mid-0.5*D],[0 1024]);

% for right branch
% subplot(132);imagesc(im);colormap(gray)
% line([mid+701 mid+701],[0 1024]);line([L+701 L+701],[0 1024]);line([R+701 R+701],[0 1024]);
% D=80;
% subplot(133);imagesc(ImMeanR);
% line([mid mid],[0 1024]);line([mid+0.5*D mid+0.5*D],[0 1024]);line([mid-0.5*D mid-0.5*D],[0 1024]);
