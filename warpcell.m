function warpcell

% When warpcell is run in a folder of aerofracture RAW .bmp files, the
% image sequence will be transformed such that the view on the bead packing
% is completely top-down

Imgs=struct(dir('*.bmp'));
[List{1:length(Imgs)}]=deal(Imgs.name);
if isempty(Imgs)
    Imgs=struct(dir('*.png'));
    [List{1:length(Imgs)}]=deal(Imgs.name);
end
mkdir('Projected');
%I=imread(List{1});
m=msgbox('Import the RAW initial image for this sequence');
uiwait(m);
data=uiimport('-file');
image=fields(data);
I=getfield(data,image{:});
h=figure;imshow(I);
Wcell=32;% packing width (cm)
Lcell=inputdlg('What is the length of the packing in cm?', 'Length',1);
Lcell=str2num(Lcell{1}); % packing length (cm)
R=imrect;r=getPosition(R);
%R=createMask(R);
uiwait(h);
r=[r(1),r(2);r(1)+r(3),r(2);r(1)+r(3),r(2)+r(4);r(1),r(2)+r(4)];
h=figure;imshow(I);
ROI=impoly;roi=getPosition(ROI);
ROI=createMask(ROI);
uiwait(h);
tform = fitgeotrans(roi,r, 'projective');
%P=imwarp(ROI.*double(I),tform)/255;
ROI2=imwarp(ROI,tform);
[ry,rx]=meshgrid(1:size(ROI2,1),1:size(ROI2,2));
rx=rx';ry=ry';
xmin=min(min(rx(ROI2(:)>0)));ymin=min(min(ry(ROI2(:)>0)));
width=max(max(rx(ROI2(:)>0)))-xmin;height=max(max(ry(ROI2(:)>0)))-ymin;
%P=imcrop(P,[xmin,ymin,width,height]);
Ap=height*width;
Am=Wcell*Lcell;
pixpercm=sqrt(Ap/Am); % pixels per cm
%
h=waitbar(0,'Projecting images');
for i=1:length(Imgs)
    I=(imread(List{i}));
    P=imwarp(ROI.*double(I),tform)/255;
    P=im2uint8(P);
    P=imcrop(P,[xmin,ymin,width,height]);
    imwrite(P,['Projected/Image_' sprintf('%05d',i) '.png'],'png')
    waitbar(i/length(Imgs),h,'Projecting images');
end
waitbar(1,h,'Done')
close(h);
save('Projected/dims.mat','pixpercm','Wcell','Lcell')