%% Plot the interaction time distribution in log-log scale
% load the sorted bee or human data

dtlist = round(durations); 
tlog = logspace(log10(min(dtlist))-0.1,log10(max(dtlist))+0.1,30)'; % constant logarithmic spacing of the plot
% plus and minus 0.1 is to include the minimum and maximum interaction time
% in the range 

tlogbinN = size(tlog,1)-1; % number of bins
tcounts = zeros([tlogbinN 1]); % the counts of the bin
terr = zeros([tlogbinN 1]); % standard error of the mean of each bin
tlogx = zeros([tlogbinN 1]);
for i=1:tlogbinN
    tpart = dtlist(tlog(i)<=dtlist); % partition the list of interaction times according to the logarithmic spacing
    tpart2 = tpart(tpart<tlog(i+1));
    tcounts(i) = size(tpart2,1)/(tlog(i+1)-tlog(i)); % logarithmic binning
    tlogx(i) = (tlog(i+1)+tlog(i))/2; % take the middle value of the bin for the x-value
    terr(i) = sqrt(size(tpart2,1))/(tlog(i+1)-tlog(i));
end

figure()
loglog(tlogx,tcounts./sum(tcounts),'.','MarkerSize',12)
hold on
te1 = errorbar(tlogx,tcounts./sum(tcounts),terr./sum(tcounts)); % normalized
set(get(te1,'Parent'),'xscale','log','yscale','log')
te1.LineStyle = 'None';
te1.LineWidth = 1.5;
xlabel('t (s)')
ylabel('P(t)')
hold off


%% plot the interaction time distribution in log-linear scale
%{
dtlist = round(durations); 
tlin = linspace(min(dtlist)-0.1,max(dtlist)+0.1,30)'; % constant logarithmic spacing of the plot
tlinbinN = size(tlin,1)-1; % number of bins
tlincounts = zeros([tlinbinN 1]); % the counts of the bin
terr = zeros([tlinbinN 1]); % standard error of the mean of each bin
tlinx = zeros([tlinbinN 1]);
for i=1:tlinbinN
    tpart = dtlist(tlin(i)<=dtlist);
    tpart2 = tpart(tpart<tlin(i+1));
    tlincounts(i) = size(tpart2,1)/(tlin(i+1)-tlin(i));
    tlinx(i) = (tlin(i+1)+tlin(i))/2; % take the middle value of the bin
    terr(i) = sqrt(size(tpart2,1))/(tlin(i+1)-tlin(i));
end

figure()
semilogy(tlinx,tlincounts./sum(tlincounts),'.','MarkerSize',12)
hold on
xlabel('t (s)')
ylabel('P(t)')
hold off
%}

%% Record the pair interaction information

beepairlist = [bee1s bee2s];
pairs = unique(beepairlist,'rows'); % unique pairs of bees
numpairs = size(pairs,1); % number of pairs
numdata = size(beepairlist,1); % length of the total data (number of total interactions)

% make a cell to divide the raw data according to pairs
fptk = cell(numpairs, 2); % {n,1}: pair [i,j], {n,2}: list of pair interaction times
freq = zeros(numpairs,1); % number of interactions for each pair
mfpt = zeros(numpairs,1); % mean of pair interaction times

comp = pairs(1,1:2); % first pair to compare
numInt = 0; % number of interactions
pnum = 1; % index to indicate which pair
fptk{pnum,1} = comp;

for i=1:numdata % go through the sorted data
    if beepairlist(i,1:2)==comp % if the pair is still that pair
        numInt = numInt + 1;
        fptk{pnum,2} = [fptk{pnum,2} durations(i)]; % add duration of the interaction into the list
    else % the pair changed
        freq(pnum) = numInt; % record the number of contacts of the previous pair
        mfpt(pnum) = sum(fptk{pnum,2})/numInt; % calculate the mean interacton time for the pair
        numInt = 1; % initialize
        pnum = pnum+1; % next pair
        comp = pairs(pnum,1:2); % next pair to compare
        fptk{pnum,1} = comp;
        fptk{pnum,2} = durations(i);
    end
end
freq(pnum) = numInt; % record the number of contacts for the last pair
mfpt(pnum) = sum(fptk{pnum,2})/numInt; % mean pair interaction time of the last pair


%% check whether f(w,t) is exponential - test the pair with the most interaction times
[maxfreq,maxfi] = max(freq); % maximum freq and its index
maxfreqtlist = fptk{maxfi,2}; % interaction times of the most freq pair
maxfreqtau = mfpt(maxfi); % mean interaction time of the most freq pair

fwtpd = makedist('Exponential',maxfreqtau);

figure()
qqplot(maxfreqtlist,fwtpd) % qq-plot to test the exponentiality


%% logarithmic binning of mean time

mfpt = round(mfpt);
mlog = logspace(log10(min(mfpt))-0.1,log10(max(mfpt))+0.1,30)'; % constant logarithmic spacing of the plot
mlogbinN = size(mlog,1)-1; % number of bins
mcounts = zeros([mlogbinN 1]); % the counts of the bin
merr = zeros([mlogbinN 1]); % standard error of the mean of each bin
mlogx = zeros([mlogbinN 1]);
for i=1:mlogbinN
    mpart = mfpt(mlog(i)<=mfpt); % partition the mean interaction time according to the log spacing
    mpart2 = mpart(mpart<mlog(i+1));
    mcounts(i) = size(mpart2,1)/(mlog(i+1)-mlog(i)); % logarithmic binning
    mlogx(i) = (mlog(i+1)+mlog(i))/2; % take the middle value of the bin
    merr(i) = sqrt(size(mpart2,1))/(mlog(i+1)-mlog(i));
end

figure()
loglog(mlogx,mcounts./sum(mcounts),'.','MarkerSize',12)
hold on
me1 = errorbar(mlogx,mcounts./sum(mcounts),merr./sum(mcounts));
set(get(me1,'Parent'),'xscale','log','yscale','log')
me1.LineStyle = 'None';
me1.LineWidth = 1.5;
xlabel('\tau (s)')
ylabel('P(\tau)')
hold off
%% Plot the interaction time and mean interaction time distributions together
figure()
loglog(tlogx,tcounts./sum(tcounts),'k.','MarkerSize',12)
hold on
loglog(mlogx,mcounts./sum(mcounts),'rx','MarkerSize',8,'LineWidth',2)

te1 = errorbar(tlogx,tcounts./sum(tcounts),terr./sum(tcounts));
set(get(te1,'Parent'),'xscale','log','yscale','log')
te1.LineStyle = 'None';
te1.LineWidth = 1.5;
set(te1,'color','k')

me1 = errorbar(mlogx,mcounts./sum(mcounts),merr./sum(mcounts));
set(get(me1,'Parent'),'xscale','log','yscale','log')
me1.LineStyle = 'None';
me1.LineWidth = 1.5;
set(me1,'color','r')

legend('Interaction time','Mean pair time','location','best')
xlabel('t (s)')
ylabel('P(t)')
hold off

%% Save the processed_data, which I actually mean the analyzed data obtained above.
save('datasetName_processed_data_date.mat','bee1s','bee2s','durations','mfpt','freq','fptk','tlogx','tcounts','terr','mlogx','mcounts','merr')



%% showing the same power for t distribution and tau distribution

% Here we distinguish between tau_ij and tau_i. tau_ij is the mean pair
% interaction time associated with an energy barrier and is obtained from
% the fitting parameter. On the other hand, tau_i is the actual average of pair
% interaction times, and this is what mfpt was. 

% load processed_data and fit data

% interaction time distribution
dtlist = round(durations); 
tlog = logspace(log10(min(dtlist))-0.1,log10(max(dtlist))+0.1,30)'; % constant logarithmic spacing of the plot
tlogbinN = size(tlog,1)-1; % number of bins
tcounts = zeros([tlogbinN 1]); % the counts of the bin
tlogx = zeros([tlogbinN 1]);
terr = zeros([tlogbinN 1]);
for i=1:tlogbinN
    tpart = dtlist(tlog(i)<=dtlist);
    tpart2 = tpart(tpart<tlog(i+1));
    tcounts(i) = size(tpart2,1)/(tlog(i+1)-tlog(i));
    tlogx(i) = (tlog(i+1)+tlog(i))/2; % take the middle value of the bin
    terr(i) = sqrt(size(tpart2,1))/(tlog(i+1)-tlog(i));
end

% mean pair interaction time distribution (tau_ij)
w1_mt7 = hecoeff_w1f;
w2_mt7 = hecoeff_w2f;

% getting rid of the artifact outlier in datasets
w2_mt7s = sort(w2_mt7); 
w2_mt7 = w2_mt7s(2:end); % jha1:(2:end), jha2:(3:end), f2f:(3:end), timtroph4:(2:end), timtroph5:(2:end)

wcollect = [w1_mt7;w2_mt7]; % get the mean interaction time from fitting.
mdur = 1./wcollect; % tau is 1/w

mlog = logspace(log10(min(mdur))-0.1,log10(max(mdur))+0.1,30)'; % constant logarithmic spacing of the plot
mlogbinN = size(mlog,1)-1; % number of bins
mcounts = zeros([mlogbinN 1]); % the counts of the bin
mlogx = zeros([mlogbinN 1]);
merr = zeros([mlogbinN 1]);
for i=1:mlogbinN
    mpart = mdur(mlog(i)<=mdur); % logarithmic partition
    mpart2 = mpart(mpart<mlog(i+1));
    mcounts(i) = size(mpart2,1)/(mlog(i+1)-mlog(i)); % logarithmic binning
    mlogx(i) = (mlog(i+1)+mlog(i))/2; % take the middle value of the bin
    merr(i) = sqrt(size(mpart2,1))/(mlog(i+1)-mlog(i));
end

% mean pair interaction time (tau_i)
lmlist = round(mfpt); 
lmlog = logspace(log10(min(mfpt))-0.1,log10(max(mfpt))+0.1,30)'; % constant logarithmic spacing of the plot
lmlogbinN = size(lmlog,1)-1; % number of bins
lmcounts = zeros([lmlogbinN 1]); % the counts of the bin
lmlogx = zeros([lmlogbinN 1]);
lmerr = zeros([lmlogbinN 1]);
for i=1:lmlogbinN
    lmpart = lmlist(lmlog(i)<=lmlist); % logarithmic partition
    lmpart2 = lmpart(lmpart<lmlog(i+1));
    lmcounts(i) = size(lmpart2,1)/(lmlog(i+1)-lmlog(i)); % logarithmic binning
    lmlogx(i) = (lmlog(i+1)+lmlog(i))/2; % take the middle value of the bin
    lmerr(i) = sqrt(size(lmpart2,1))/(lmlog(i+1)-lmlog(i));
end

%% Save the data
save('datasetName_scalinginfo_date.mat','tcounts','terr','tlogx','mcounts','merr','mlogx','lmcounts','lmerr','lmlogx')

%% Plot f(t), p(tau_ij), p(tau_i) together
figure()
loglog(tlogx,tcounts./sum(tcounts),'ro','MarkerSize',5,'MarkerFaceColor','r')
hold on
loglog(mlogx,mcounts./sum(mcounts),'bx','MarkerSize',7,'LineWidth',2)
loglog(lmlogx,lmcounts./sum(lmcounts),'ks','MarkerSize',3,'MarkerFaceColor','k')

te = errorbar(tlogx,tcounts./sum(tcounts),terr./sum(tcounts));
set(get(te,'Parent'),'xscale','log','yscale','log')
te.LineStyle = 'None';
te.LineWidth = 1.5;
set(te,'color','r')

me = errorbar(mlogx,mcounts./sum(mcounts),merr./sum(mcounts));
set(get(me,'Parent'),'xscale','log','yscale','log')
me.LineStyle = 'None';
me.LineWidth = 1.5;
set(me,'color','b')

lme = errorbar(lmlogx,lmcounts./sum(lmcounts),lmerr./sum(lmcounts));
set(get(lme,'Parent'),'xscale','log','yscale','log')
lme.LineStyle = 'None';
lme.LineWidth = 1.5;
set(lme,'color','k')

box on
xlabel('');
ylabel('');
% xlim([10^(0) 10^2])
% ylim([10^(-6) 1])
set(gca,'XTick',[10^0 10^2]);
set(gca,'TickLength',[0.05 0.05])
set(gca,'YTick',[10^(-4) 10^(-2) 10^0])
set(gca,'TickLength',[0.05 0.05])
hold off



%% Histogram of the weight in the hyperexponential fit

numpf = length(hecoeff_pf); % number of fitted pairs
pf_bal = zeros(numpf,1); % this is to make sure that the weight >= 0.5.
% this constraint is there because weight = 0.2 and 0.8 would represent the
% same weight.
for j=1:numpf
    if hecoeff_pf(j)>0.5
        pf_bal(j) = hecoeff_pf(j);
    else
        pf_bal(j) = 1-hecoeff_pf(j);
    end
end

figure()
histogram(pf_bal)
hold on
xlabel('weight on one exponential')
ylabel('counts')
hold off
