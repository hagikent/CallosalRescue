%% SpontanEventProcess

clear all

save_dir='Z:\AnalysisData\hagihara\KirGCaMPMurakami\150403\';

dir1='Z:\AnalysisData\hagihara\KirGCaMPMurakami\150403\xyt001_ch1_rp';
dir2='Z:\AnalysisData\hagihara\KirGCaMPMurakami\150403\xyt002_ch1_rp';
%{
dir3='Z:\AnalysisData\hagihara\141010\xyt004_ch1_rp';
dir4='Z:\AnalysisData\hagihara\141010\xyt005_ch1_rp';
dir5='Z:\AnalysisData\hagihara\141010\xyt006_ch1_rp';

dir6='Z:\AnalysisData\hagihara\141010_2\xyt002_ch1_rp';
dir7='Z:\AnalysisData\hagihara\141010_2\xyt003_ch1_rp';
dir8='Z:\AnalysisData\hagihara\141010_2\xyt004_ch1_rp';
dir9='Z:\AnalysisData\hagihara\141010_2\xyt005_ch1_rp';
dir10='Z:\AnalysisData\hagihara\141010_2\xyt006_ch1_rp';
dir11='Z:\AnalysisData\hagihara\141010_2\xyt007_ch1_rp';
dir12='Z:\AnalysisData\hagihara\141010_2\xyt008_ch1_rp';

dir13='Z:\AnalysisData\hagihara\141010_3\xyt001_ch1_rp';
dir14='Z:\AnalysisData\hagihara\141010_3\xyt002_ch1_rp';
dir15='Z:\AnalysisData\hagihara\141010_3\xyt003_ch1_rp';
dir16='Z:\AnalysisData\hagihara\141010_3\xyt004_ch1_rp';
dir17='Z:\AnalysisData\hagihara\141010_3\xyt005_ch1_rp';
dir18='Z:\AnalysisData\hagihara\141010_3\xyt006_ch1_rp';
dir19='Z:\AnalysisData\hagihara\141010_3\xyt007_ch1_rp';
dir20='Z:\AnalysisData\hagihara\141010_3\xyt008_ch1_rp';

dir21='Z:\AnalysisData\hagihara\141011\xyt002_ch1_rp';
dir22='Z:\AnalysisData\hagihara\141011\xyt003_ch1_rp';
dir23='Z:\AnalysisData\hagihara\141011\xyt004_ch1_rp';
dir24='Z:\AnalysisData\hagihara\141011\xyt005_ch1_rp';
dir25='Z:\AnalysisData\hagihara\141011\xyt006_ch1_rp';
dir26='Z:\AnalysisData\hagihara\141011\xyt007_ch1_rp';
dir27='Z:\AnalysisData\hagihara\141011\xyt008_ch1_rp';
dir28='Z:\AnalysisData\hagihara\141011\xyt009_ch1_rp';
dir29='Z:\AnalysisData\hagihara\141011\xyt010_ch1_rp';

dir30='Z:\AnalysisData\hagihara\141011_2\xyt001_ch1_rp';
dir31='Z:\AnalysisData\hagihara\141011_2\xyt002_ch1_rp';
dir32='Z:\AnalysisData\hagihara\141011_2\xyt003_ch1_rp';
dir33='Z:\AnalysisData\hagihara\141011_2\xyt004_ch1_rp';
dir34='Z:\AnalysisData\hagihara\141011_2\xyt005_ch1_rp';
dir35='Z:\AnalysisData\hagihara\141011_2\xyt006_ch1_rp';
dir36='Z:\AnalysisData\hagihara\141011_2\xyt007_ch1_rp';
%{
dir37='M:\AnalysisData\hagihara\130601_3\Image6_ch2_rp';
dir38='M:\AnalysisData\hagihara\130601_3\Image7_ch2_rp';
dir39='M:\AnalysisData\hagihara\130601_3\Image9_ch2_rp';
dir40='M:\AnalysisData\hagihara\130601_3\Image10_ch2_rp';
dir41='M:\AnalysisData\hagihara\130601_3\Image12_ch2_rp';

dir42='M:\AnalysisData\hagihara\130601_4\scan001_ch1_rp';
dir43='M:\AnalysisData\hagihara\130601_4\scan002_ch1_rp';
dir44='M:\AnalysisData\hagihara\130601_4\scan003_ch1_rp';
dir45='M:\AnalysisData\hagihara\130601_4\scan004_ch1_rp';
dir46='M:\AnalysisData\hagihara\130601_4\scan005_ch1_rp';
%}
%}
animal1=[1:2];
%animal2=[6:12];
%animal3=[13:20];
%animal4=[21:29];
%animal5=[30:36];

runN=2;
AnimalN=1;
Lcriteria=0.20;
HLcriteria=0.60;

animal0=0;   %do not touch
flag.save=1;
flag.PosiNega=1;
%%
animal0=0;   %do not touch
Eventpartrate=[];
Cellpartrate=[];
Posi2Hrate=[];
Posi2Lrate=[];
Posi2Krate=[];
Nega2Hrate=[];
Nega2Lrate=[];
Nega2Krate=[];
HFreq=[];
LFreq=[];
KFreq=[];

PosiNum=[];

if flag.save == 1
    mkdir(save_dir)
end


for ii=1:runN
    dir=eval(['dir',num2str(ii)]);
    cd(dir)
    load eventmatrix2
    load Recmin
    if flag.PosiNega==1
    load PosiNum
    end
    Rectime(ii)=Recmin;
    celln(ii)=size(eventmatrix2,1);
    eventn(ii)=size(eventmatrix2,2);
    NegaNum=[1:celln(ii)];
    NegaNum(PosiNum)=[];    
    EPR=sum(eventmatrix2,1)./celln(ii);
    CPR=sum(eventmatrix2,2)./eventn(ii);
    
    
    Heventnum=find(EPR>HLcriteria);
    Keventnum=find(EPR<Lcriteria);
    Leventnum=[1:eventn(ii)];
    Leventnum([Heventnum,Keventnum])=[];
    
    Hfreq=length(Heventnum)./Recmin;
    Lfreq=length(Leventnum)./Recmin;
    Kfreq=length(Keventnum)./Recmin;
    
    Hnum(ii)=length(Heventnum);
    Lnum(ii)=length(Leventnum);
    Knum(ii)=length(Keventnum);
    
    Hevent=eventmatrix2(:,Heventnum);
    Levent=eventmatrix2(:,Leventnum);
    Kevent=eventmatrix2(:,Keventnum);
    
    P2H=sum(Hevent(PosiNum,:),2)./size(Hevent,2);
    P2L=sum(Levent(PosiNum,:),2)./size(Levent,2);
    P2K=sum(Kevent(PosiNum,:),2)./size(Kevent,2);
    N2H=sum(Hevent(NegaNum,:),2)./size(Hevent,2);
    N2L=sum(Levent(NegaNum,:),2)./size(Levent,2);
    N2K=sum(Kevent(NegaNum,:),2)./size(Kevent,2);
    
    P2Hfreq{ii}=sum(Hevent(PosiNum,:),2)./Recmin;
    P2Lfreq{ii}=sum(Levent(PosiNum,:),2)./Recmin;
    P2Kfreq{ii}=sum(Kevent(PosiNum,:),2)./Recmin;
    N2Hfreq{ii}=sum(Hevent(NegaNum,:),2)./Recmin;
    N2Lfreq{ii}=sum(Levent(NegaNum,:),2)./Recmin;
    N2Kfreq{ii}=sum(Kevent(NegaNum,:),2)./Recmin;
    
    Eventpartrate=[Eventpartrate EPR];
    Cellpartrate=[Cellpartrate ;CPR];
    Posi2Hrate=[Posi2Hrate ;P2H];
    Posi2Lrate=[Posi2Lrate ;P2L];
    Posi2Krate=[Posi2Krate ;P2K];
    Nega2Hrate=[Nega2Hrate ;N2H];
    Nega2Lrate=[Nega2Lrate ;N2L];
    Nega2Krate=[Nega2Krate ;N2K];
    HFreq=[HFreq ;Hfreq];
    LFreq=[LFreq ;Lfreq];
    KFreq=[KFreq ;Kfreq];
end

for ii=1:runN
    if Hnum(ii)==0
         P2Hfreq{ii}=[];
         N2Hfreq{ii}=[];
    end
    if Lnum(ii)==0
         P2Lfreq{ii}=[];
         N2Lfreq{ii}=[];
    end
    if Knum(ii)==0
         P2Kfreq{ii}=[];
         N2Kfreq{ii}=[];
    end
end


%%

for ii=1:AnimalN
    animalind=eval(['animal',num2str(ii)]);
if sum(HFreq(animalind))==0;
    HFreqAnimalAve(ii)=0;
    HFreqAnimalStd(ii)=0;
else
    HFreqAnimalAve(ii)=sum(HFreq(animalind))./length(animalind);
    HFreqAnimalStd(ii)=std(HFreq(animalind))./sqrt(length(animalind));
end
    LFreqAnimalAve(ii)=sum(LFreq(animalind))./length(animalind);
    LFreqAnimalStd(ii)=std(LFreq(animalind))./sqrt(length(animalind));
    KFreqAnimalAve(ii)=sum(KFreq(animalind))./length(animalind);
    KFreqAnimalStd(ii)=std(KFreq(animalind))./sqrt(length(animalind));
end

%%
eventnummatrix(1,:)=Hnum;
eventnummatrix(2,:)=Lnum;
eventnummatrix(3,:)=Knum;

%% plot
timesum=sum(Rectime)
bins10=histc([Eventpartrate],[0:0.1:1])./timesum;

figure(1),bar(bins10);


if flag.save==1;
    cd(save_dir)
    saveas(figure(1),'bins10hist','fig');
end
%%
figure(2),
subplot(2,3,1),hist(Posi2Hrate),title('Kir Posi to H event')
subplot(2,3,2),hist(Posi2Lrate),title('Kir Posi to L event')
subplot(2,3,3),hist(Posi2Krate),title('Kir Posi to K event')
subplot(2,3,4),hist(Nega2Hrate),title('Kir Nega to H event')
subplot(2,3,5),hist(Nega2Lrate),title('Kir Nega to L event')
subplot(2,3,6),hist(Nega2Krate),title('Kir Nega to K event')

if flag.save==1;
    saveas(figure(2),'PosiNega2KLMevent','fig')
end
%%
figure(3),
subplot(1,2,1),
plot(cumsum(histc(Posi2Hrate,[0:0.1:1]))*100/max(cumsum(histc(Posi2Hrate,[0:0.1:1.1]))),'r','LineWidth',3),hold on
plot(cumsum(histc(Nega2Hrate,[0:0.1:1]))*100/max(cumsum(histc(Nega2Hrate,[0:0.1:1]))),'LineWidth',3)
xlim([0,11])
subplot(1,2,2),
plot(cumsum(histc(Posi2Lrate,[0:0.1:1]))*100/max(cumsum(histc(Posi2Lrate,[0:0.1:1.1]))),'r','LineWidth',3),hold on
plot(cumsum(histc(Nega2Lrate,[0:0.1:1]))*100/max(cumsum(histc(Nega2Lrate,[0:0.1:1]))),'LineWidth',3)
xlim([0,11])
%%
figure(4),
barweb([HFreqAnimalAve ;LFreqAnimalAve ;KFreqAnimalAve],[HFreqAnimalStd ;LFreqAnimalStd ;KFreqAnimalStd])
ylim([0,3.75])
tcksX={'H';'L';'K';};set(gca, 'XTickLabel',tcksX ,'XTick',1:length(tcksX));
if flag.save==1;
    saveas(figure(4),'HLKeventFreq','fig')
end

%%
P2HfreqAll=cell(AnimalN,1);
P2LfreqAll=cell(AnimalN,1);
P2KfreqAll=cell(AnimalN,1);
N2HfreqAll=cell(AnimalN,1);
N2LfreqAll=cell(AnimalN,1);
N2KfreqAll=cell(AnimalN,1);

%%
for jj=1:AnimalN
for ii=max(eval(['animal' num2str(jj-1)]))+1:max(eval(['animal' num2str(jj)]))
    foldernum=ii-max(eval(['animal' num2str(jj-1)]));
    P2HfreqAll{jj}=[P2HfreqAll{jj} ;P2Hfreq{ii}];
    P2LfreqAll{jj}=[P2LfreqAll{jj} ;P2Lfreq{ii}];
    P2KfreqAll{jj}=[P2KfreqAll{jj} ;P2Kfreq{ii}];
    N2HfreqAll{jj}=[N2HfreqAll{jj} ;N2Hfreq{ii}];
    N2LfreqAll{jj}=[N2LfreqAll{jj} ;N2Lfreq{ii}];
    N2KfreqAll{jj}=[N2KfreqAll{jj} ;N2Kfreq{ii}];
end
end
%%
figure(5),
subplot(2,3,1),bar(histc([eval('P2HfreqAll{:}')],[0:0.2:1.2]),1,'r'),title('Kir Posi to H event'),tcksX={'-0.2';'-0.4';'-0.6';'-0.8';'-1.0';'-1.2';'1.2'};set(gca, 'XTickLabel',tcksX ,'XTick',1:length(tcksX));
subplot(2,3,2),bar(histc([eval('P2LfreqAll{:}')],[0:0.2:1.2]),1,'r'),title('Kir Posi to L event'),set(gca, 'XTickLabel',tcksX ,'XTick',1:length(tcksX));
subplot(2,3,3),bar(histc([eval('P2KfreqAll{:}')],[0:0.2:1.2]),1,'r'),title('Kir Posi to K event'),set(gca, 'XTickLabel',tcksX ,'XTick',1:length(tcksX));
subplot(2,3,4),bar(histc([eval('N2HfreqAll{:}')],[0:0.2:1.2]),1),title('Kir Nega to H event'),set(gca, 'XTickLabel',tcksX ,'XTick',1:length(tcksX));
subplot(2,3,5),bar(histc([eval('N2LfreqAll{:}')],[0:0.2:1.2]),1),title('Kir Nega to L event'),set(gca, 'XTickLabel',tcksX ,'XTick',1:length(tcksX));
subplot(2,3,6),bar(histc([eval('N2KfreqAll{:}')],[0:0.2:1.2]),1),title('Kir Nega to K event'),set(gca, 'XTickLabel',tcksX ,'XTick',1:length(tcksX));

if flag.save==1;
    saveas(figure(5),'freq of cells','fig')
end




%%

figure(6),
subplot(1,2,1),
plot(cumsum(histc([eval('P2HfreqAll{:}')],[0:0.2:1.2]))*100/max(cumsum(histc([eval('P2HfreqAll{:}')],[0:0.1:1.2]))),'r','LineWidth',2),hold on
plot(cumsum(histc([eval('N2HfreqAll{:}')],[0:0.2:1.2]))*100/max(cumsum(histc([eval('N2HfreqAll{:}')],[0:0.1:1.2]))),'LineWidth',2)
set(gca, 'XTickLabel',tcksX ,'XTick',1:length(tcksX));
subplot(1,2,2),
plot(cumsum(histc([eval('P2LfreqAll{:}')],[0:0.2:1.2]))*100/max(cumsum(histc([eval('P2LfreqAll{:}')],[0:0.1:1.2]))),'r','LineWidth',2),hold on
plot(cumsum(histc([eval('N2LfreqAll{:}')],[0:0.2:1.2]))*100/max(cumsum(histc([eval('N2LfreqAll{:}')],[0:0.1:1.2]))),'LineWidth',2)
set(gca, 'XTickLabel',tcksX ,'XTick',1:length(tcksX));
%%
mean(Posi2Lrate(find(Posi2Lrate>-1)));
mean(Nega2Lrate(find(Nega2Lrate>-1)));
mean(Posi2Hrate(find(Posi2Hrate>-1)));
mean(Nega2Hrate(find(Nega2Hrate>-1)));

median(Posi2Lrate(find(Posi2Lrate>-1)));
median(Nega2Lrate(find(Nega2Lrate>-1)));
median(Posi2Hrate(find(Posi2Hrate>-1)));
median(Nega2Hrate(find(Nega2Hrate>-1)));
%%
if flag.save==1;
    save freqAll P2HfreqAll P2LfreqAll P2KfreqAll N2HfreqAll N2LfreqAll N2KfreqAll
end

%%
if flag.PosiNega==1
    
[H_Leve,P_eveL]=kstest2(Posi2Lrate(find(Posi2Lrate>-1)),Nega2Lrate(find(Nega2Lrate>-1)))
[H_Heve,P_eveH]=kstest2(Posi2Hrate(find(Posi2Hrate>-1)),Nega2Hrate(find(Nega2Hrate>-1)))

PLrs=ranksum(Posi2Lrate(find(Posi2Lrate>-1)),Nega2Lrate(find(Nega2Lrate>-1)))
PHrs=ranksum(Posi2Hrate(find(Posi2Hrate>-1)),Nega2Hrate(find(Nega2Hrate>-1)))

end
%%

