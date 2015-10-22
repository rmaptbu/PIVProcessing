function[mid,w]=FindMidlineCalib(prof)
% profs=abs(smooth(prof,5));
% plot(1:length(prof),-prof,'.-');hold all

% find approximate positions of edges
[sizesl,edgesl]=findpeaks(abs(double((prof(1:275)))),'MINPEAKHEIGHT',500,'MINPEAKDISTANCE',3);   % find left edge
[~,Il]=sort(sizesl);      % put in order of size
ledge=edgesl(Il(end));    % select two largest peaks

% find approximate positions of edges
[sizesr,edgesr]=findpeaks(double(abs(prof(276:end))),'MINPEAKHEIGHT',500,'MINPEAKDISTANCE',3);   % find left edge
[~,Ir]=sort(sizesr);      % put in order of size
redge=edgesr(Ir(end))+275;    % select two largest peaks

Buffer=20;
if length(sizesl)==1
    [~,L]=max(abs(prof(ledge-Buffer:ledge+Buffer)));L=L+ledge-1-Buffer;
elseif length(sizesl)>1
    a=edgesl(Il(end-1));b=edgesl(Il(end));
    if a<b A=a;B=b;else A=b;B=a;end
    [~,L]=min(prof(A:B));L=L+A-1;
end

if length(sizesr)==1
    [~,R]=max(abs(prof(redge-Buffer:redge+Buffer)));R=R+redge-1-Buffer;
elseif length(sizesr)>1
    a=edgesr(Ir(end-1));b=edgesr(Ir(end));
    if a<b A=a;B=b;else A=b;B=a;end
    [~,R]=max(prof((A:B)+275));R=R+275+A-1;
end

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
