% Analysis of Spontaneous Activity (in developing cortex)
%flags plot:plot wave of all cells
%      celltypesep:separate posi and nega (F+/F-)
%      
%120508 Hagihara.K
%130819 modified for rawtimecourses event process 

clear all
close all
flag.plot=0; % be careful, very slow...
flag.cellrecog=0; % 0:do cellrecog/1:load labelimg (already exist)
flag.CellTypeSep=1;  % 0:no-sep/1:separate manually /2:load Posinum(already exist)
flag.NEPmethod=1; % the way of deNoising contamination from neuropil  1:NEP / 0:BV (recommend 1)
DN_switch=1;  %  0;remove baseline when calculating DNcoef/1;don't care (when flag.NEPmethod=0)
testplot=1;  % 
%% set dir&file&rec time
anal_dir='Z:\AnalysisData\hagihara\141010\xyt002_ch1_rp_test';              %set dir
data_dir='Z:\RawData\hagihara\Nikon\141010\mat\';      %set dir
 
template=16;
cell_size=12;
%template=20;
%cell_size=18;

Recmin=10;        %set duration       
Recsec=Recmin*60;

filename1='xyt002_ch1_rp'; 

fn=1;
%%
for iiii=1:fn

    close all
    
    keep flag DN_switch testplot  anal_dir data_dir  filename  template Recmin Recsec...
        filename1 filename2 filename3 filename4 filename5 filename6 filename7Å@filename8 filename9...
        fn iiii
    
    filename=eval(['filename',num2str(iiii)]);

%%
cd(anal_dir)
mkdir(cd,filename)
cd(filename)
dir_save=[anal_dir,'\',filename];

load([data_dir,filename]); %xyt data is loaded as stack
avgimg_fname=[data_dir,filename,'_avg.tif'];
img=readtiff(avgimg_fname); %xy for cell recog (average of xyt data)

array=stack;
clear stack
[Y,X,T]=size(array);   % pixel*pixel*frame


fps=T./Recsec;
save Recmin Recmin
save fps fps
save frameNum T

%% param
if flag.cellrecog==0
%params for cell recognition
%%template=16;
r_th=0.5;
%cell_size=12;

save params4cellrecog template r_th cell_size
end
%% All ROI timecourse
Avetimecourse=squeeze(mean(mean(array,1),2));
Avetimecourse_norm=Avetimecourse./mean(Avetimecourse);
Avetimecourse_norm=tcLowCut(Avetimecourse_norm,size(Avetimecourse_norm,1)/4,'gaussian',1);
figure(1);plot(Avetimecourse_norm);hold on
%plot([1:size(Avetimecourse_norm)],2*std(Avetimecourse_norm)+1,'-k');
hold off

saveas(figure(1),'Avetimecourse_norm','fig')
save Avetimecourse Avetimecourse
save Avetimecourse_norm Avetimecourse_norm

%% cell recog
if flag.cellrecog==0
[labelimg,hp_img]=imFindCellsTM(img,template, r_th,cell_size, 1);

%% 
save labelimg labelimg
save hp_img hp_img

im=[imShade(hp_img, logical(labelimg)),repmat(hp_img-0.5,[1,1,3])];
figure,
imshow(im);
imwrite(im,'shadecells.tif');

%%
if flag.NEPmethod==0
%% get BV image

disp('Get BV image!!')

a_th_BV2_low=15;
a_th_BV2_high=1000;
intns_th=0.7;

input_BV=1;

while input_BV==1

a_th_BV2_low0=a_th_BV2_low;
a_th_BV2_low=input(['>>Pleae type new params; a_th_low (current value: ', num2str(a_th_BV2_low) ,', default)>>' ]) ;
    if isempty(a_th_BV2_low)
        a_th_BV2_low=a_th_BV2_low0;
        clear a_th_BV2_low0
    else
        clear a_th_BV2_low0
    end

a_th_BV2_high0=a_th_BV2_high;
a_th_BV2_high=input(['>>Pleae type new params; a_th_high (current value: ', num2str(a_th_BV2_high) ,', default)>>' ]) ;
    if isempty(a_th_BV2_high)
        a_th_BV2_high=a_th_BV2_high0;
        clear a_th_BV2_high0
    else
        clear a_th_BV2_high0
    end
    
intns_th0=intns_th;
intns_th=input(['>>Pleae type new params; intns_th (current value: ', num2str(intns_th) ,', default)>>' ]) ;
    if isempty(intns_th)
        intns_th=intns_th0;
        clear intns_thr0
    else
        clear intns_thr0
    end
            
hp_bv=hp_img<intns_th;

labelimg_BV2=bwlabel(hp_bv);
tempstats=regionprops(labelimg_BV2, 'Area');
A=[tempstats.Area];
ind=find(A< a_th_BV2_low);
    if ~isempty(ind)
            for e=1:length(ind)
                labelimg_BV2(labelimg_BV2==ind(e))=labelimg_BV2(labelimg_BV2==ind(e))*0;
            end
    %        labelimg=bwlabel(labelimg); 
    %        clear ind 
    end

ind=find(A> a_th_BV2_high);
    if ~isempty(ind)
            for e=1:length(ind)
                labelimg_BV2(labelimg_BV2==ind(e))=labelimg_BV2(labelimg_BV2==ind(e))*0;
            end
    %        labelimg=bwlabel(labelimg); 
    %        clear ind 
    end

labelimg_BV=sort_labelimg(labelimg_BV2); 

im_BV=[imshade(hp_img*-1, logical(labelimg_BV)), imshade(hp_img, logical(labelimg_BV))];
figure(300);hold on;clf(300);imshow(im_BV);

in=input('>>Use this labeling? Yes: type 1, otherwise, set new params (default). >>');
    if isempty(in)
        in = 0; 
    else
    end
    
    if in==1
        input_BV=0;
    else
    end
    
end    

save params4BVrecog a_th_BV2_low  a_th_BV2_high  intns_th
%% %%%%%%%%% manual removal
n=0;
while n==0
        figure(gcf); hold on;
        im=[imShade_ty(hp_img*-1, zeros(size(labelimg_BV)),logical(labelimg_BV)),imShade(hp_img-0.5, logical(labelimg_BV)),repmat(hp_img-0.5,[1,1,3])];
        imshow(im);
        
        disp('>>Please click labeled BVs you want to remove, and press Enter.>>')
        [x y]=ginput;  % select position to be removed
       labelimg_manualRm=labelimg_BV;

       if ~isempty(x)  
                
                Ncells_RM=zeros(size(x,1),3);
                ii=0;
                     for j=1:size(x,1)
                         tempX=round(x(j));
                         tempY=round(y(j));
                         tempN=labelimg_manualRm(tempY, tempX);
                            if tempN==0
                                 disp('>>Cell was not found!');
                            else
                             labelimg_manualRm(labelimg_manualRm==tempN)=labelimg_manualRm(labelimg_manualRm==tempN)*0;
                             im_BV=[imShade_ty(hp_img*-1, zeros(size(labelimg_manualRm)),logical(labelimg_manualRm)),imShade(hp_img-0.5, logical(labelimg)),repmat(hp_img-0.5,[1,1,3])];
                             figure(gcf);clf(gcf);hold on;imshow(im_BV);
                             ii=ii+1;
                             Ncells_RM(ii,1)=tempN;
                             Ncells_RM(ii,2)=tempX;
                             Ncells_RM(ii,3)=tempY;
                            end
                     end
                     Ncells_RM=Ncells_RM(1:ii,:);
                
            else
                Ncells_RM=[];
            end
             in2=input('>>Finish manual removal or Retry? Finish: 1, Retry: 0 (default 0)>>');
                if isempty(in2)
                     in2=0;
                 end

                 if in2==1
                     labelimg_BV=sort_labelimg(labelimg_manualRm);
                     clear labelimg_manualRm
                     
%                      cd(fullfile(analysis_dir,'labelimg'));
                      cd(dir_save)
                      save('labelimg_BV', 'labelimg_BV', 'hp_img');
                      save cellrecogParams_BV a_th_BV2_low a_th_BV2_high 
                      imwrite(im_BV,'labelimg_BV.tif');
                    input_BV=0;
                     
                     
                     n=1;
                     
                 clear labelimg_manualRm Ncells_RM
                 im=[imShade(hp_img*-1, logical(labelimg_BV)),imShade(hp_img-0.5, logical(labelimg_BV)),repmat(hp_img-0.5,[1,1,3])];
                 figure(gcf);clf;hold on;imshow(im);
                 end
     
                
end  

tempstatsBV=regionprops(labelimg_BV,'EquivDiameter');
r_BVs=round([tempstatsBV.EquivDiameter]);
%%
end
end

if flag.cellrecog==1
    load labelimg
    %load params4cellrecog
    load hp_img
end
%% compute OGB timecourses

timeCourses = stackGetTimeCourses(array, labelimg);  
save timeCourses timeCourses

celln=size(timeCourses,2);
for ii=1:celln
    timeCourses_norm(:,ii)=timeCourses(:,ii)./mean(timeCourses(:,ii));
end

save timeCourses_norm timeCourses_norm


%% compute BV timecourse
if flag.NEPmethod==0;
timeCourses_BV = stackGetTimeCourses(array, labelimg_BV);

save timeCourses_BV timeCourses_BV
end
%% make neuropil time courses around cells
rs_factor=1;
r_cells=ceil(template/2);
ringSize=round(rs_factor*r_cells)+1;
excludeSize=2; 
minPixels=5;
% intensity_th=1.05;

[NeuropilMask, NeuropilLabel] = getNeuropilMask(labelimg, ringSize, excludeSize, minPixels);
im_np=[imShade_ty(hp_img, logical(labelimg), logical(NeuropilLabel)),imShade_ty(hp_img, (zeros(size(labelimg))), logical(NeuropilLabel)),repmat(hp_img-0.5,[1,1,3])];
figure;imshow(im_np)

cd(dir_save)
save NPlabelparams ringSize excludeSize minPixels
save NeuropilLabel NeuropilLabel
writetiff8(NeuropilLabel.*10000000, 'NeuropilLabel.tif', 1);
writetiff8(im_np, 'Cell_NeuropilLabel.tif', 1);

NeuropilTimeCourses = getNeuropilTimeCourses(array, NeuropilMask);
cd(dir_save)
save NeuropilTimeCourses NeuropilTimeCourses



%% NP timeCourse around BV
if flag.NEPmethod==0;

ringSize_BVNP=round(max(round(r_BVs)))*rs_factor+1;
[BVMask, BVLabel] = getNeuropilMask(labelimg_BV, ringSize-1, excludeSize, minPixels);

NPLabelBV=logical(BVLabel);
im_bv=[imShade_ty(hp_img*-1, logical(labelimg_BV ), logical(NPLabelBV)),imShade_ty(hp_img, (logical(labelimg)), logical(NPLabelBV)),repmat(hp_img-0.5,[1,1,3])];

cd(dir_save)
figure;imshow(im_bv)
saveas(gcf, 'NPlabelonBV')

cd(dir_save)
save NPlabelBVparams ringSize_BVNP excludeSize minPixels
save BVLabel BVLabel
writetiff8(BVLabel.*10000000, 'BVLabel.tif', 1);
writetiff8(im_bv, 'BV_NeuropilLabel.tif', 1);

BVMaskmod=BVMask;
if length(size(BVMask))==2
    BVMaskmod(:,:,1)=BVMaskmod;
    BVMaskmod(:,:,2)=BVMaskmod;
    timeCourses_BV(:,2)=timeCourses_BV;
end

NeuropilTimeCourses_BV = getNeuropilTimeCourses(array, BVMaskmod);
cd(dir_save)
saveas(gcf, 'NPlabelonBV')

%% compute denoised coefficient (BV method)

NCcoefBV=zeros(size(NeuropilTimeCourses_BV,2),1);
cd(dir_save)

for i=1:size(NeuropilTimeCourses_BV,2);
 if DN_switch==1   
    NCcoefBV(i)=mean(timeCourses_BV(:,i),1)/mean(NeuropilTimeCourses_BV(:,i),1);
 else
    NCcoefBV(i)=(NeuropilTimeCourses_BV(:,i)-mean(NeuropilTimeCourses_BV(:,i)))\(timeCourses_BV(:,i)-mean(timeCourses_BV(:,i)));
 end
end

coefBV=mean(NCcoefBV(NCcoefBV>0))

cd(dir_save)
save coefBV coefBV NCcoefBV
    
for i=1:size(NeuropilTimeCourses,2);
if DN_switch==0
deNeuropilTimeCourses(:,i)=timeCourses(:,1)-((NeuropilTimeCourses(:,i)-mean(NeuropilTimeCourses(:,i))).*coefBV);
elseif DN_switch==1
deNeuropilTimeCourses(:,i)=timeCourses(:,i)-NeuropilTimeCourses(:,i).*coefBV; 
end
deNeuropilTimeCourses_norm(:,i)=deNeuropilTimeCourses(:,i)./mean(deNeuropilTimeCourses(:,i));
NeuropilTimeCourses_norm(:,i)=NeuropilTimeCourses(:,i)./mean(NeuropilTimeCourses(:,i));
end

cd(dir_save)
save deNeuropilTimeCourses deNeuropilTimeCourses
save deNeuropilTimeCourses_norm deNeuropilTimeCourses_norm
save NeuropilTimeCourses_norm NeuropilTimeCourses_norm
end
%% Neuropil SignalÅ@DeNoising (Non Event Period:NEP Method)
n=10;
Wn=1.0/fps;  %params for high cut (tcHighCut)    ex) fps=3.65(Galvo)Å®Wn=0.4

if fps<1
    Wn=0.8;
end

for ii=1:celln
   
nep=find(timeCourses(:,ii)<mean(timeCourses(:,ii)));

Ymogeraw=timeCourses(nep,ii);
Ymoge=tcHighCut(Ymogeraw,n,Wn);
Ymoge=tcLowCut(Ymoge,size(Ymoge,1)/2,'gaussian',1);
Ymoge=Ymoge-mean(Ymoge);

Xmogeraw=NeuropilTimeCourses(nep,ii);
Xmoge=tcHighCut(Xmogeraw,n,Wn);
Xmoge=tcLowCut(Xmoge,size(Xmoge,1)/2,'gaussian',1);
Xmoge=Xmoge-mean(Xmoge);

Xvec=[sum(Xmoge.^2),sum(Xmoge);sum(Xmoge),length(Xmoge)];
Yvec=[sum(Xmoge.*Ymoge);sum(Ymoge)];
Avec(:,ii)=Xvec\Yvec;                 % y=Ax+B  A=Avec(1),B=Avec(2)=0,least squares

end

%%
for ii=1:celln
    NPTC=NeuropilTimeCourses(:,ii)-mean(NeuropilTimeCourses(:,ii));
    deNeuropilTimeCourses2(:,ii)=timeCourses(:,ii)-NPTC.*Avec(1,ii)-Avec(2,ii);
    deNeuropilTimeCourses2_norm(:,ii)=deNeuropilTimeCourses2(:,ii)./mean(deNeuropilTimeCourses2(:,ii));
end

save deNeuropilTimeCourses2 deNeuropilTimeCourses2
save deNeuropilTimeCourses2_norm deNeuropilTimeCourses2_norm
%% Test Plot
if testplot==1
%%    
plotmoge=10;
figure,
subplot(4,1,1);plot(timeCourses(:,plotmoge));title('timeCourse')
subplot(4,1,2);plot(NeuropilTimeCourses(:,plotmoge));title('Neuropil')
subplot(4,1,3);plot((NeuropilTimeCourses(:,plotmoge)-mean(NeuropilTimeCourses(:,plotmoge))).*Avec(1,plotmoge)+mean(NeuropilTimeCourses(:,plotmoge)));
subplot(4,1,4);plot(deNeuropilTimeCourses2(:,plotmoge));title('deNoised')

%%    
end


%% LCÅ®LP
if flag.NEPmethod==1
deNeuropilTimeCourses=deNeuropilTimeCourses2;
end

%LC
LCcutoff=5; %min
LCTC=tcLowCut(deNeuropilTimeCourses,size(deNeuropilTimeCourses,1)/(Recmin/LCcutoff),'gaussian',1);

%LP
for ii=1:celln
LPLCTC(:,ii)=tcHighCut(LCTC(:,ii),n,Wn);
LPLCTC_norm(:,ii)=LPLCTC(:,ii)./mean(LPLCTC(:,ii));
end


%% LCÅ®LP for raw time courses
%LC
LCcutoff=5; %min
LCTCraw=tcLowCut(timeCourses,size(timeCourses,1)/(Recmin/LCcutoff),'gaussian',1);

%LP
for ii=1:celln
LPLCTCraw(:,ii)=tcHighCut(LCTCraw(:,ii),n,Wn);
LPLCTCraw_norm(:,ii)=LPLCTCraw(:,ii)./mean(LPLCTCraw(:,ii));
end



%% calculate base line noise(SD) emulating Konnerth lab's Anal.

%BaseLineNoise=deNeuropilTimeCourses_norm-smoothedTimeCourses1_norm;
%SD=std(BaseLineNoise);
%save SD SD

%BaseLineNoise2=timeCourses_norm-smoothedTimeCourses0_norm;
%SD2=std(BaseLineNoise2);
%save SD2 SD2
%% plotting LPÅ®LC
%Criteria=1.1
%plotn=17;
%figure,
%subplot(3,1,1),plot(timeCourses(:,plotn)),hold on
%               plot(NeuropilTimeCourses(:,plotn).*coefBV,'r');
%               plot(deNeuropilTimeCourses(:,plotn),'g'); hold off
%               title([num2str(plotn),'     Blue:raw Red:Neuropil*contamination coef G:denoisedTC'])
              
%subplot(3,1,2),plot(LPTC(:,plotn));hold on
%               plot(deNeuropilTimeCourses(:,plotn),'g');hold off
%               title('LowPassFiltered')
            
%subplot(3,1,3),plot(LCLPTC_norm(:,plotn)); hold on
%               plot([1:size(LCLPTC_norm,1)],Criteria)
%               title('LowCutFiltered')
%               ylim([0.8,1.6]);

%% plotting2  
Criteria=1.05
save Criteria Criteria
plotn=20;

switch flag.NEPmethod
    case 0
figure,
subplot(3,1,1),plot(timeCourses(:,plotn)),hold on
               plot(NeuropilTimeCourses(:,plotn).*coefBV,'r');
               plot(deNeuropilTimeCourses(:,plotn),'g'); hold off
               title([num2str(plotn),'     Blue:raw Red:Neuropil*contamination coef G:denoisedTC'])
              
subplot(3,1,2),plot(LCTC(:,plotn));hold on
               plot(deNeuropilTimeCourses(:,plotn),'g');hold off
               title('LowCutFiltered')
            
subplot(3,1,3),plot(LPLCTC_norm(:,plotn)); hold on
               plot([1:size(LPLCTC_norm,1)],Criteria)
               title('LowPassFiltered')
               ylim([0.8,1.6]);               
    case 1
figure,
subplot(3,1,1),plot(timeCourses(:,plotn)),hold on
               plot((NeuropilTimeCourses(:,plotn)-mean(NeuropilTimeCourses(:,plotn))).*Avec(1,plotn),'r');
               plot(deNeuropilTimeCourses(:,plotn),'g'); hold off
               title([num2str(plotn),'     Blue:raw Red:Neuropil*contamination coef G:denoisedTC'])
              
subplot(3,1,2),
               plot(deNeuropilTimeCourses(:,plotn),'g');hold on
               plot(LCTC(:,plotn));hold off
               title('LowCutFiltered')
            
subplot(3,1,3),plot(LPLCTC_norm(:,plotn)); hold on
               plot([1:size(LPLCTC_norm,1)],Criteria)
               title('LowPassFiltered')
               ylim([0.8,1.2]);  
end

%% raster plot

raster1=zeros(celln,T);

for ii=1:celln
    for jj=1:T
        if LPLCTC_norm(jj,ii)>Criteria
            raster1(ii,jj)=1;
        else
            raster1(ii,jj)=0;
        end
    end
end

%% raster plot for raw

raster1raw=zeros(celln,T);

for ii=1:celln
    for jj=1:T
        if LPLCTCraw_norm(jj,ii)>Criteria
            raster1raw(ii,jj)=1;
        else
            raster1raw(ii,jj)=0;
        end
    end
end

%% event process 
convwindow=[1 1 1]; % window for convolution. if you want to fill in N size gaps,use N+1 size window 

raster=raster1;

eventsum=sum((raster),1);
eventsumbw=logical(eventsum);
eventsumbw=conv(eventsumbw,convwindow);   %% filling in the gaps. Size of the gaps are defined by the convwindow size. 
eventsumbw=eventsumbw(1:T);
eventsumbw=logical(eventsumbw);
eventlabel=bwlabel(eventsumbw);

eventmatrix=zeros(celln,max(eventlabel));

for ii=1:max(eventlabel)
    for jj=1:celln
        if sum(raster(jj,find(eventlabel==ii)))>0
            eventmatrix(jj,ii)=1;
        else
            eventmatrix(jj,ii)=0;
        end
    end
end
%figure(31),subplot(2,1,2),imshow(imcomplement(eventsumbw))
figure(32),imshow(imcomplement(eventmatrix));
save eventmatrix eventmatrix
saveas(figure(32),'eventmatrix','fig')


%% event process for raw
convwindow=[1 1 1]; % window for convolution. if you want to fill in N size gaps,use N+1 size window 

raster=raster1raw;

eventsumraw=sum((raster),1);
eventsumbw=logical(eventsumraw);
eventsumbw=conv(eventsumbw,convwindow);   %% filling in the gaps. Size of the gaps are defined by the convwindow size. 
eventsumbw=eventsumbw(1:T);
eventsumbw=logical(eventsumbw);
eventlabel=bwlabel(eventsumbw);

eventmatrixraw=zeros(celln,max(eventlabel));

for ii=1:max(eventlabel)
    for jj=1:celln
        if sum(raster(jj,find(eventlabel==ii)))>0
            eventmatrixraw(jj,ii)=1;
        else
            eventmatrixraw(jj,ii)=0;
        end
    end
end
%figure(31),subplot(2,1,2),imshow(imcomplement(eventsumbw))
figure(64),imshow(imcomplement(eventmatrixraw));
save eventmatrixraw eventmatrixraw
saveas(figure(64),'eventmatrixraw','fig')

%% denoise raster and eventmatrix (>20% participation rate)

eventper=eventsum./celln;
raster2=raster;
for ii=1:T
    if eventper(ii)<0.05
            raster2(:,ii)=0;
    end
end

event2sum=sum((raster2),1);
event2sumbw=logical(event2sum);
event2sumbw=conv(event2sumbw,convwindow);   %% filling in the gaps. Size of the gaps are defined by the convolution window size. 
event2sumbw=event2sumbw(1:T);
event2sumbw=logical(event2sumbw);
event2label=bwlabel(event2sumbw);

eventmatrix2=zeros(celln,max(event2label));

for ii=1:max(event2label)
    for jj=1:celln
        if sum(raster2(jj,find(event2label==ii)))>0
            eventmatrix2(jj,ii)=1;
        else
            eventmatrix2(jj,ii)=0;
        end
    end
end

figure(33),subplot(2,1,1),imshow(imcomplement(raster2))
subplot(2,1,2),imshow(imcomplement(event2sumbw))
set(gca,'Position',get(gca,'OuterPosition'));
figure(34),imshow(imcomplement(eventmatrix2))
set(gca,'Position',get(gca,'OuterPosition'));
save event2label event2label
save raster2 raster2
save eventmatrix2 eventmatrix2
saveas(figure(33),'raster2','fig')
saveas(figure(34),'eventmatrix2','fig')

%% denoise raster and eventmatrix (>20% participation rate) for raw

eventperraw=eventsumraw./celln;
raster2raw=raster1raw;
for ii=1:T
    if eventperraw(ii)<0.20
            raster2raw(:,ii)=0;
    end
end

event2sum=sum((raster2raw),1);
event2sumbw=logical(event2sum);
event2sumbw=conv(event2sumbw,convwindow);   %% filling in the gaps. Size of the gaps are defined by the convolution window size. 
event2sumbw=event2sumbw(1:T);
event2sumbw=logical(event2sumbw);
event2labelraw=bwlabel(event2sumbw);

eventmatrix2raw=zeros(celln,max(event2labelraw));

for ii=1:max(event2labelraw)
    for jj=1:celln
        if sum(raster2raw(jj,find(event2labelraw==ii)))>0
            eventmatrix2raw(jj,ii)=1;
        else
            eventmatrix2raw(jj,ii)=0;
        end
    end
end

figure(66),subplot(2,1,1),imshow(imcomplement(raster2raw))
subplot(2,1,2),imshow(imcomplement(event2sumbw))
set(gca,'Position',get(gca,'OuterPosition'));
figure(68),imshow(imcomplement(eventmatrix2raw))
set(gca,'Position',get(gca,'OuterPosition'));
save event2labelraw event2labelraw
save raster2raw raster2raw
save eventmatrix2raw eventmatrix2raw
saveas(figure(66),'raster2raw','fig')
saveas(figure(68),'eventmatrix2raw','fig')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%cell type separate (manual)%%%%%%%%
if flag.CellTypeSep~=0;
%%
switch flag.CellTypeSep
    case 1
PosiNum=getCellNumbers('labelimg', avgimg_fname);
      
PosiNum=sort(PosiNum);
saveas(gcf,'selectedcells.fig');

PosiNum= PosiNum(PosiNum>0);
save PosiNum PosiNum
    case 2
        load('Posinum')
end
        
%%

eventmatrixT=eventmatrix;
eventmatrixT(PosiNum,:)=0;

eventmatrixColor(:,:,1)=imcomplement(eventmatrixT);
eventmatrixColor(:,:,2)=imcomplement(eventmatrix);
eventmatrixColor(:,:,3)=imcomplement(eventmatrixT);
eventmatrixColor(PosiNum,:,:)=eventmatrixColor(PosiNum,:,:).*0.8;
figure(35),imshow(eventmatrixColor);
set(gca,'Position',get(gca,'OuterPosition'));
saveas(figure(35),'eventmatrixColor','fig')

eventmatrix2T=eventmatrix2;
eventmatrix2T(PosiNum,:)=0;

eventmatrix2Color(:,:,1)=imcomplement(eventmatrix2T);
eventmatrix2Color(:,:,2)=imcomplement(eventmatrix2);
eventmatrix2Color(:,:,3)=imcomplement(eventmatrix2T);
eventmatrix2Color(PosiNum,:,:)=eventmatrix2Color(PosiNum,:,:).*0.8;
figure(36),imshow(eventmatrix2Color);
saveas(figure(36),'eventmatrix2Color','fig')

raster1posi=raster;
raster1posi(PosiNum,:)=0;

raster1Color(:,:,1)=imcomplement(raster1posi);
raster1Color(:,:,2)=imcomplement(raster);
raster1Color(:,:,3)=imcomplement(raster1posi);

raster1Color(PosiNum,:,:)=raster1Color(PosiNum,:,:).*0.8;

figure(37),imshow(raster1Color);
saveas(figure(37),'raster1Color','fig')

raster2posi=raster2;
raster2posi(PosiNum,:)=0;

raster2Color(:,:,1)=imcomplement(raster2posi);
raster2Color(:,:,2)=imcomplement(raster2);
raster2Color(:,:,3)=imcomplement(raster2posi);

raster2Color(PosiNum,:,:)=raster2Color(PosiNum,:,:).*0.8;
div(:,:,1)=zeros(10,T);
div(:,:,2)=zeros(10,T);
div(:,:,3)=ones(10,T).*0.5;
mark(:,:,1)=imcomplement([event2sumbw;event2sumbw;event2sumbw;event2sumbw;event2sumbw;event2sumbw;event2sumbw;event2sumbw;event2sumbw]);
mark(:,:,2)=imcomplement([event2sumbw;event2sumbw;event2sumbw;event2sumbw;event2sumbw;event2sumbw;event2sumbw;event2sumbw;event2sumbw]);
mark(:,:,3)=imcomplement([event2sumbw;event2sumbw;event2sumbw;event2sumbw;event2sumbw;event2sumbw;event2sumbw;event2sumbw;event2sumbw]);

raster2Color=[raster2Color;div;mark];

figure(38),imshow(raster2Color);
saveas(figure(38),'raster2Color','fig')
%%
NegaNum=[1:size(raster,1)];
NegaNum(PosiNum)=[];

rasterColorPNAS(:,:,1)=[zeros(size(raster2(PosiNum,:))); raster2(NegaNum,:)];
rasterColorPNAS(:,:,2)=[raster2(PosiNum,:); raster2(NegaNum,:)];
rasterColorPNAS(:,:,3)=[zeros(size(raster2(PosiNum,:))); raster2(NegaNum,:)];
save rasterColorPNAS rasterColorPNAS
%%
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot all raw wave %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag.plot==1
%%   
ScreenSize=get(0,'ScreenSize');
figure('Position',[0 0 ScreenSize(3).*0.98 ScreenSize(4).*0.90])
switch exist('PosiNum')
    case 1
for ii=1:celln
    if sum(ii==PosiNum)==1
        color=[1 0 0];
    else
        color=[0 0 1];
    end
    subplot(ceil(celln./10),10,ii); 
    plot(LPLCTC_norm(:,ii),'Color',color); hold on
    %plot([1:T],1.,'r');
    plot([1:T],1.05,'g');
    title(num2str(ii));
     ylim([0.8,1.3]);  
     xlim([0,T]);
end
    case 0
for ii=1:celln
    subplot(ceil(celln./10),10,ii); 
    plot(LPLCTC_norm(:,ii)); hold on
    %plot([1:T],1.1,'r');
    plot([1:T],1.05,'g');
    title(num2str(ii));
     ylim([0.8,1.3]);  
     xlim([0,T]);
end        
end
%% %%%%%%%%%%%%%%%%%%%%%%%  semitra  %%%%%%%%%%%%%%%%%%%%%%%%
ScreenSize=get(0,'ScreenSize');
figure('Position',[0 0 ScreenSize(3).*0.98 ScreenSize(4).*0.90]),
for ii=1:celln
    xx=[1:T]; yy=LPLCTC_norm(xx,ii);
    p(ii)=pplot(xx,yy); hold on
    set(p(ii),'EdgeAlpha',0.05);
    set(p(ii),'EdgeColor',[rand(1) rand(1) rand(1)]);
    set(p(ii),'LineWidth',0.5);
end
    plot([1:T],1.1,'r');
    plot([1:T],1.2,'g');hold off

%% posi plot 
figure,
for ii=1:size(PosiNum,1)
    subplot(size(PosiNum,1),1,ii);
    plot(timeCourses_norm(:,PosiNum(ii))); hold on
    %plot(NeuropilTimeCourses_norm(:,PosiNum(ii)).*coefBV,'r');
    plot(LPLCTC_norm(:,PosiNum(ii)),'g');hold off
    title(['PosiCellsTC       ' num2str(PosiNum(ii)), '   B:Raw       G:deNoise+filter'])
end

%%
figure,
for ii=1:size(PosiNum,1)
    subplot(size(PosiNum,1),1,ii);
    plot(timeCourses(:,PosiNum(ii))); hold on
    %plot(NeuropilTimeCourses(:,PosiNum(ii)).*coefBV,'r');
    plot(LPLCTC(:,PosiNum(ii)),'g');hold off
    title(['PosiCellsTC       ' num2str(PosiNum(ii)), '      B:Raw        G:deNoise+filter'])
end

%%
end
end
        



