function aerobox2(BW,foldername)

% aerobox(BW) performs box counting on the binary aerofracture image BW.
% A .mat file is saved with the resulting box counts and size information.

%Test of empty image
if sum(BW(:))<1
    Db=0;
    cnt=NaN;
    sid=NaN;
    save('aerobox.mat','Db','sid','cnt');
    return
end
%Estimate largest box dimension
[X,Y]=meshgrid(1:size(BW,2),1:size(BW,1));
X=fliplr(X);
lf=max(max(double(BW>0).*X));
if (lf/2)~=round(lf/2)
    lf=lf+1;
end
yl=max(Y(:));
if (yl/2)~=round(yl/2)
    yl=yl+1;
end
if lf>=yl
    pad=ceil([(lf/2)-(yl/2),0]);
    BW2=padarray(BW,pad,0,'both');
    x0=size(BW2,2)-lf;
    y0=1;
    sidecheck=1;
    BW2=BW2(y0:size(BW2,1),x0:size(BW2,2));
    BW2=padarray(BW2,[max(size(BW2))-size(BW2,1) max(size(BW2))-size(BW2,2)],0,'post');
    lf=size(BW2,1);
else
    pad=ceil([0,(yl/2)-(lf/2)]);
    BW2=padarray(BW,pad,0,'both');
    x0=size(BW2,2)-yl;
    y0=1;
    sidecheck=0;
    BW2=BW2(y0:y0+yl-1,x0:x0+yl-1);
    BW2=BW2(y0:size(BW2,1),x0:size(BW2,2));
    BW2=padarray(BW2,[max(size(BW2))-size(BW2,1) max(size(BW2))-size(BW2,2)],0,'post');
    yl=size(BW2,1);
end

%Find lower cutoff
m=msgbox(['Import the "Local_frac_vars" mat-file from ' foldername]);
uiwait(m);
data=uiimport('-file');
Af=data.Af;
X=data.X;
Af=Af{length(Af)}; %if cell
m=msgbox('Find typical finger thickness from plot, input in command window');
uiwait(m);
plot(X,Af,'b.');
if sidecheck==1
    lowcut=round(2*(lf/input('Waiting for typical finger thickness : ')));
else
    lowcut=round(2*(yl/input('Waiting for typical finger thickness : ')));
    if lowcut<round(2*(yl/floor(lf/4)))
        lowcut=round(2*(yl/floor(lf/4)));
    end
end
pause(0.1);
close all

%Perform box counting
h=waitbar(0,'Counting.. please wait');
%n=1;
%idx=1;
for n=1:lowcut
    %while n<=lf
    %cnt(n)=0;
    if sidecheck==1
        sid(n)=(lf/n);
        [x,y]=meshgrid(1:n,0:n:n*(n-1));
        indexed_boxes=imresize((x+y),[lf lf],'nearest');
        cnt(n)=length(unique(indexed_boxes.*double(BW2)))-1;
        waitbar(n/lowcut,h,'Counting.. please wait');
        %n=n+(5^(floor(log10(n)))); %Size dependent logarithmic number step
        %idx=idx+1;
    else
        sid(n)=(yl/n);
        [x,y]=meshgrid(1:n,0:n:n*(n-1));
        indexed_boxes=imresize((x+y),[yl yl],'nearest');
        cnt(n)=length(unique(indexed_boxes.*double(BW2)))-1;
        waitbar(n/lowcut,h,'Counting.. please wait');
        %n=n+(5^(floor(log10(n)))); %Size dependent logarithmic number step
        %idx=idx+1;
    end
    
end
waitbar(1,h,'Done');
close(h);
sidall=sid;cntall=cnt;
plot(log10(sidall)
sid=sid(cnt~=((1:lowcut).^2)); %remove lower cutoff values
cnt=cnt(cnt~=((1:lowcut).^2)); %remove lower cutoff values
% When yl>lf
if sidecheck==0
    cnt=cnt(sid<(lf/4));
    sid=sid(sid<(lf/4));
end
[~,Db]=linreg_coeffs(log10(sid),log10(cnt)); %Find slope
Db=-round(Db*100)/100;
save('aerobox.mat','sid','cnt','Db','sidall','cntall');

