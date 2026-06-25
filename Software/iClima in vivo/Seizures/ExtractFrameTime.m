%%Extract the frame time of the imaging form the electrophysiology

%Fai partire la prima parte del programma e otterrai un cell array che
%contiene tutte le finestre di tempo di acquisizione del 2P, comprese le
%partenze errate. puoi osservare nel plot quelle che vuoi eliminare ed
%immetterle nell'array discardsTw. IMPORTANTE, per ottenere il corretto
%numero di frame l'app prende il segnale dalla banda passante 1 che
%dovrebbe essere settata tra 0 e 60Hz in modo da rendere il segnale più
%smooth altrimenti il numero di volte che la treshold viene superata
%potrebbe variare per via dell'oscillazione nel segnale del ch2.
app=zebra %inserire come handle zebra
a=app.bandPassed_LFP(1,2,1,:);
a=squeeze(a)';
b=zebra.workLFP(1,1,:);
b=squeeze(b)';
xt=app.x;
%plot(app.x,a)
tresholdWindow=-10; %always check with the plot if makes sense
tresholdFrame=0;
frameLog=a>tresholdFrame;
frameDiff=diff(frameLog);
frameDiff=[0,frameDiff];
windowLog=a>tresholdWindow;
windowDiff=diff(windowLog);
windowDiff=[0,windowDiff];
frameIndex=find(frameDiff==-1);
frameTime=xt(frameIndex);
windowIndxUp=find(windowDiff==1);
windowIndxDown=find(windowDiff==-1);
for i=1:length(windowIndxUp)
    timeWindowsIndx{i}= [windowIndxUp(i) windowIndxDown(i)];
    VisualizationTimeWindows{i}=xt( windowIndxDown(i))-xt(windowIndxUp(i));%not used in the code, just to see better which are the window to eliminate
end

FrameCutTime=[];
%Now check the timewindow and decide which one you want to discard
discardsTw=[ 10 11 12 15 16];
timeWindowsIndx(discardsTw) = [];
indxArray=cat(2, timeWindowsIndx{:});

for i=1:2:length(indxArray)
    FrameCutTime=[FrameCutTime, frameTime(frameTime>=xt(indxArray(i)) & frameTime<=xt(indxArray(i+1)))];
end
    

%Guarda qui per vedere se ci sono frame sbagliati, in caso diminuire la
%banda passante.
figure
plot(diff([0 FrameCutTime]),'o-') %ti plotta la differenza nel tempo, quindi se due frame sono troppo vicini lo mostra.



% writematrix(FrameCutTime','FrameCutTime.txt');
% writematrix(b','LFP.txt');
%writematrix(xt', 'xt.txt');


