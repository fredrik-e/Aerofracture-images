function longfrac

Imgs=struct(dir('*.png'));
[List{1:length(Imgs)}]=deal(Imgs.name);
i=1;
I=im2bw(imread(List{i}));
m=size(I,1);
n=size(I,2);
[X,Y]=meshgrid(1:n,1:m);
X=fliplr(X);
h=waitbar(0,'Estimating pattern properties');
for i=1:length(Imgs)
    I=im2bw(imread(List{i}));
    longfinger(i)=max(max(I.*X));
    mass(i)=sum(I(:));
    ytmp=I.*Y;
    if max(ytmp(:))>0
        Yrange(i)=max(ytmp(:))-min(ytmp(ytmp>0));
    else
        Yrange(i)=0;
    end
    waitbar(i/length(Imgs),h,'Estimating pattern properties');
end
waitbar(1,h,'Done');
close(h);
save('Area-lengths.mat','longfinger','mass','Yrange','m','n');