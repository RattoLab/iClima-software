%% Processing of the imaging of seizure experiment
%  Allineare immagini 
%1) prendi un frame di tamplate per la funzione "myregistration"
%2) fai un average (frame by frame) di verde e rosso
%3) registri con "myregistration" la tua stack, con il tuo template
%4) prendi i deltas di myregistration e le passi "imreg_trasl" aprendo
%prima sono il canale verde e poi solo il canale rosso da registrare. In
%questo modo ti assicuri anche che i due canali siano allineati entrambi
%allo stesso modo

%Load the Ch1 channel
%Nota che il nome delle immagini verde e rossa devono essere le stesse
[file_in,tmpdir_in] = uigetfile();
info = imfinfo([tmpdir_in,file_in]);
app.n_wl = numel(info);
app.y= info(1).Height;
app.x= info(1).Width;
app.gre=nan(app.y,app.x,app.n_wl);
for i=1:app.n_wl
    app.gre(:,:,i) = imread([tmpdir_in,file_in],"tif",i);
end

% import ch2
file_in_red= strrep(file_in,'Ch1','Ch2');
app.red=nan(app.y,app.x,app.n_wl);
for i=1:app.n_wl
    app.red(:,:,i) = imread([tmpdir_in,file_in_red],"tif",i);
end

%Average of Green and Red channel 
for i=1:app.n_wl
    tmp=cat(3,app.gre(:,:,i),app.red(:,:,i));
    app.avg(:,:,i) = mean(tmp,3);
end
%take deltas
[registered,deltas]= myregistration(app.avg,app.avg(:,:,1));
%register green channel
registered_Ch1 = imreg_trasl(app.gre,deltas);
%register Red channel
registered_Ch2 = imreg_trasl(app.red,deltas);
%prepare to save the image
            tagstruct.Photometric=Tiff.Photometric.MinIsBlack;
            tagstruct.Compression=Tiff.Compression.PackBits;
            tagstruct.BitsPerSample=32;
            tagstruct.SamplesPerPixel=1;
            tagstruct.SubFileType=2;
            tagstruct.SampleFormat=Tiff.SampleFormat.IEEEFP;
            tagstruct.ImageLength=app.y;
            tagstruct.ImageWidth=app.x;
            tagstruct.PlanarConfiguration=Tiff.PlanarConfiguration.Chunky;
            tagstruct.TileLength=app.y;
            tagstruct.TileWidth=app.x;
%save Ch1
t = Tiff([tmpdir_in 'registered_Ch1.tif'],'w');
setTag(t,tagstruct)
t.write(single(registered_Ch1(:,:,1)))
for i=2:app.n_wl
    writeDirectory(t);
    setTag(t,tagstruct)
    t.write(single(registered_Ch1(:,:,i)));
end
close(t)

%save Ch2
t = Tiff([tmpdir_in 'registered_Ch2.tif'],'w');
setTag(t,tagstruct)
t.write(single(registered_Ch2(:,:,1)))
for i=2:app.n_wl
    writeDirectory(t);
    setTag(t,tagstruct)
    t.write(single(registered_Ch2(:,:,i)));
end
close(t)

%% Take values from eROIClo
%carica manualmente le due variabili buffer.gre e buffer.red
buffer.gre=[];
buffer.red=[];
mean.gre=[];
mean.red=[];
% max.gre=[];
% max.red=[];
thresh_mean.gre=[];
thresh_mean.red=[];
% 1.1) Correction for the relative gain change of the PMTs
g_over_r_reference=0.618;
g_over_r_measured=0.6326;
[mean.gre, mean.red]=compensate_relative_gain(buffer.gre,buffer.red,g_over_r_reference,g_over_r_measured);

%Remove Bleedtrhough
beta=0.11;      %value of beta @880nm for iClima: è il verde nel rosso
alfa=0.236;      %value of alfa @880nm for iClima: è il rosso nel verde
for i=1:size(mean.gre,1)
    for k=1:size(mean.gre,2)
        mean.gre(i,k)=(mean.gre(i,k)/(1-alfa*beta))-((mean.red(i,k)*alfa)/(1-alfa*beta));
        mean.red(i,k)= mean.red(i,k)-(beta*((mean.gre(i,k)-(alfa*mean.red(i,k)))/(1-alfa*beta)));
    end
end

%scrematura
%1) VALORE FLUORESCENZA TROPPO BASSO: setta una threshold
thresh= 100;
bin_matrix=mean.red(:,:)<thresh;
thresh_mean.gre=mean.gre;
thresh_mean.red=mean.red;
thresh_mean.gre(bin_matrix)=NaN;
thresh_mean.red(bin_matrix)=NaN;

%2) TOGLI LE ROI CHE ESCONO DAL PIANO FOCALE
%scegli la percentuale di fluorescenza della cellula sotto il quale, se la
%cellula scende, viene eliminata la sua analisi da quel frame
%più alta è la threshold più sei restrittivo con le scelte ed elimini più
%roi
percent=30;
maxima=maxk(mean.red,3);
maxima=min(maxima);
for i=1:size(mean.red,1)
    for k=1:size(mean.red,2)
        if (mean.red(i,k)*100/maxima(k))<percent
            thresh_mean.gre(i,k)=NaN;      %è la verde scremata per i due threshold
            thresh_mean.red(i,k)=NaN;      %è la rossa scremata per i due threshold
        end
    end
end
count=sum(~isnan(thresh_mean.red),2);



%process the remaining data
G_over_R= thresh_mean.gre./thresh_mean.red;
MedianInTime.gre= median(thresh_mean.gre,2,'omitnan');
MedianInTime.red= median(thresh_mean.red,2,'omitnan');
MedianInTime.GoverR= median(G_over_R,2,'omitnan');

figure
plot(MedianInTime.GoverR,'k','LineWidth', 1.5)
figure
plot(MedianInTime.gre,'g')
hold on 
plot(MedianInTime.red,'r')


norm_Gre= MedianInTime.gre/mean(MedianInTime.gre);
norm_Red= MedianInTime.red/mean(MedianInTime.red);
figure
plot(norm_Gre,'g','LineWidth', 1)
hold on 
plot(norm_Red,'r','LineWidth', 1)

%median su ROI a diversi tempi
MedianRoi.gre=median(thresh_mean.gre,1,'omitnan');
MedianRoi.red= median(thresh_mean.red,1,'omitnan');
MedianRoi.GoverR= median(G_over_R,1,'omitnan');

figure
plot(MedianRoi.GoverR,'k','LineWidth', 1.5)
figure
plot(MedianRoi.gre,'g')
hold on 
plot(MedianRoi.red,'r')

%Visualize R over R of different ROI
figure
tiledlayout(6,6)
for i=1:size(G_over_R,1)
    nexttile 
    plot(G_over_R(i,:))
    title(['ROI ' num2str(i)]) 
end
%%
figure
for i=1:size(thresh_mean.gre,2)
plot(thresh_mean.gre(:,i),'g')
hold on 
plot(thresh_mean.red(:,i),'r')
end

%%
% SALVATI thresh_mean.gre  & thresh_mean.red  & G_over_R
save('C:\Users\giak7\OneDrive\Documenti\LAB\New Chloride Sensor\mouse 2\variables','mean','thresh_mean','G_over_R','MedianInTime','count','alfa', 'beta','thresh','percent')