function [Nf,Tf,Af,X,COUNTS,LENGTHS]=locdim_and_fingerthickness(BW,linewidth)

% This function estimates the average finger thickness vs. radius as well
% as the local box counts per radius for the aerofracture in a linear cell
% where BW is the binary image

if ~exist('linewidth')
    linewidth=1;
end

xlength=size(BW,2);
ylength=size(BW,1);
yrng=[1;ylength];
n=9000;
L=(1:n)';
T=find((n./L)==round((n./L)));        %Find denominators giving equal pieces
T=T(T>0);                             %Remove 0's
S=cat(2,flipud(T),T);                 %Number of pieces & denominators
% h=waitbar(0,'Estimating, please wait..');
for i=1:xlength-linewidth+1;
    xrng=[i;i+linewidth-1];
    ptmp=improfile(BW,xrng,yrng,n,'nearest');
    %P{(xlength+1)-i}=improfile(BW,xrng,yrng,n,'nearest');
    Nf((xlength+1)-i)=sum(diff(ptmp)>0);
    Tf((xlength+1)-i)=(sum(ptmp>0)/n)*ylength;
    for j=1:size(S,1)
        counts(j)=sum(sum(reshape(ptmp,S(j,:)),1)>0,2);
        lengths(j)=ylength/S(j,2);
    end
    COUNTS{(xlength+1)-i}=counts;
    LENGTHS{(xlength+1)-i}=lengths;
%     waitbar(i/xlength,h,'Estimating, please wait..');
end
Af=Tf./Nf; %Average finger thickness
X=(1:xlength);
%waitbar(1,h,'Done!');
%close(h);
%save('Local_frac_vars.mat','Nf','Tf','Af','X','COUNTS','LENGTHS')