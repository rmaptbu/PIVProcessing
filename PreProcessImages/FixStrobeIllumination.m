function[Ga]=FixStrobeIllumination(ImMeanR,ImMask)
OverallMean=mean(ImMeanR(:))*ones(size(ImMeanR));%*1.15;
ImMeanRA=(1-ImMask).*ImMeanR+ImMask.*OverallMean; % apply mask

h=fspecial('gaussian',[25 25],3);               % smooth image with gaussian filter
ImMeanF=imfilter(ImMeanRA,h,'replicate');
% imagesc(ImMeanF);axis 

[h,w]=size(ImMeanRA);
[X,Y]=meshgrid(1:w,1:h);
X=X(:);Y=Y(:);Z=ImMeanRA(:);
% 2D gaussian fit object
gauss2=fittype(@(a1,sigmax,sigmay,x0,y0,x,y) a1*exp(-(x-x0).^2/(2*sigmax^2)-(y-y0).^2/(2*sigmay^2)),...
        'independent', {'x', 'y'},'dependent', 'z' );
a1=max(ImMeanRA(:)); % height, determine from image. may want to subtract background
sigmax=2; % guess width
sigmay=2; % guess width
x0=w/2; % guess position (center seems a good place to start)
y0=h/2; 
  % compute fit
sf=fit([X,Y],double(Z),gauss2,'StartPoint',[a1, sigmax, sigmay, x0,y0]);
G=feval(sf,[X,Y]);
[X2,Y2] = meshgrid(1:w,1:h);
Ga=zeros(size(X2));
for i=1:size(X)
    x=X(i);
    y=Y(i);
    Ga(y,x)=G(i);
end
% pcolor(X2,-Y2,ImMean);shading flat;axis image
% figure();pcolor(X2,-Y2,Ga);shading flat;axis image
Ga=Ga/max(Ga(:));

       
