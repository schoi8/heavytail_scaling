%% fit to hyperexponential distribution for a single pair with most interactions as a test
%{
[~,maxfi] = max(freq);
maxfreqtlist = fptk{maxfi,2};
[maxfreqecdf,xval]=ecdf(maxfreqtlist);
maxfreqecdf(maxfreqecdf==0)=[];
xval = xval(2:end);
%{
figure()
cdfplot(maxfreqtlist)
hold on
plot(xval,maxfreqecdf)
hold off
%}
fthe = fittype('1-a*exp(-b*x)-(1-a)*exp(-c*x)');
heop = fitoptions(fthe);
heop.Lower = [0.05 0 0];
heop.Upper = [0.95 Inf Inf];
heop.StartPoint = [0.9 0.1 0.01];
[hefit,hegof]=fit(xval,maxfreqecdf,fthe,heop);
yhe = feval(hefit,xval);

ftexp = fittype('1-exp(-a*x)');
expop = fitoptions(ftexp);
expop.Lower = 0;
expop.StartPoint = 0.1;
[expfit,expgof]=fit(xval,maxfreqecdf,ftexp,expop);
yexp = feval(expfit,xval);

figure()
plot(xval,maxfreqecdf,'.','MarkerSize',12)
hold on
plot(xval,yhe,'LineWidth',2)
plot(xval,yexp,'LineWidth',2)
%plot(xval,1-0.9*exp(-0.24*xval)-0.1*exp(-0.04*xval),'LineWidth',2)
legend('data','hyperexp','exp','location','best')
xlabel('t')
ylabel('F(t)')
hold off
%}


%% fit hyperexponential for all pairs
numpair = size(fptk,1); % number of pairs
r2pair = zeros(numpair,1); % will record R^2 values for each pair fit
toofew = 0; % just keeping track of the number of pairs
hecoeff_p = zeros(numpair,1); % weight on one of exponential terms
hecoeff_w1 = zeros(numpair,1); % rate of one of exponential terms
hecoeff_w2 = zeros(numpair,1); % rate of the other exponential term
hebee1 = zeros(numpair,1); % record whose interaction I am analyzing
hebee2 = zeros(numpair,1); 

for i=1:numpair
    pairint = fptk{i,2}; % list of pair interaction times
    [pairecdf,xval]=ecdf(pairint); % empirical cumulative distribution function
    pairecdf(pairecdf==0)=[];
    xval = xval(2:end);
    
    if length(pairecdf)>7 % minimum sample size for fitting is 3 -> 7 for better fitting
        fthe = fittype('1-p*exp(-w1*x)-(1-p)*exp(-w2*x)'); % CDF for hyperexponential distribution
        heop = fitoptions(fthe);
        maxint = max(pairint);
        heop.Upper = [0.99 Inf Inf]; % weight < 1
        heop.Lower = [0.01 0 0]; % w1 > 0, w2 > 0
        heop.StartPoint = [0.9 0.1 0.01];
        [hefit,hegof]=fit(xval,pairecdf,fthe,heop);
        
        hecoeff = coeffvalues(hefit);
        hecoeff_p(i)=hecoeff(1);
        hecoeff_w1(i)=hecoeff(2);
        hecoeff_w2(i)=hecoeff(3);
        
        beepair = fptk{i,1};
        hebee1(i)=beepair(1);
        hebee2(i)=beepair(2);
    
        r2pair(i) = hegof.rsquare;
    else
        toofew = toofew+1;
    end
end

nonzeroidx = find(r2pair); % index for non-zero elements
r2pair_f = r2pair(nonzeroidx);
hecoeff_pf = hecoeff_p(nonzeroidx);
hecoeff_w1f = hecoeff_w1(nonzeroidx);
hecoeff_w2f = hecoeff_w2(nonzeroidx);
hebee1f = hebee1(nonzeroidx);
hebee2f = hebee2(nonzeroidx);

figure()
histogram(r2pair_f)
hold on
xlabel('R^{2}')
ylabel('Counts')
hold off
%% save fitting data
save('datasetName_hefit_mt7_v1_date.mat','hebee1f','hebee2f','hecoeff_pf','hecoeff_w1f','hecoeff_w2f','hefit','hegof','r2pair_f')


%% Collapse of all data - pairfit_info
% ecdf information for each pair
numpair = size(fptk,1);
enoughint = length(hebee1f);
enoughind = 1;
pairfit_info = cell(enoughint,3); % 1: pair, 2: ecdf, 3: xval
toofew2 = 0;

for i=1:numpair
    pairint = fptk{i,2};
    [pairecdf,xval]=ecdf(pairint);
    pairecdf(pairecdf==0)=[];
    xval = xval(2:end);
    
    if length(pairecdf)>7 % 2->7
        pairfit_info{enoughind,1} = fptk{i,1};
        pairfit_info{enoughind,2} = pairecdf;
        pairfit_info{enoughind,3} = xval;
        enoughind = enoughind+1;
    else
        toofew2 = toofew2+1;
    end
end

%% Check the data collapse for the dataset
% check how good the fitting is. Plotting the fitting function value
% against ecdf. If fitting was good, there should be a cloud of data points
% around y=x reference line.

enoughleng = length(r2pair_f);
plotted = 0;

figure()
hold on
for i=1:enoughleng
    yval = pairfit_info{i,2};
    xval = pairfit_info{i,3};    
    plot(1-hecoeff_pf(i)*exp(-hecoeff_w1f(i)*xval)-(1-hecoeff_pf(i)).*exp(-hecoeff_w2f(i)*xval),yval,'.','markersize',12)
    plotted = plotted+1;
end
plot([0,1],[0,1],'k--','LineWidth',2)
xlabel('1-pe^{-w_{1}t}-(1-p)e^{-w_{2}t}')
ylabel('ECDF(t)')
hold off


%% collapse plot in terms of e^z
% for better data collapse, we defined z=w1x and rearranged the CDF so that
% only one variable is left in x-axis. 

mt7_numpair = size(pairfit_info,1); % number of pairs used for fitting
mt7_expz = cell(mt7_numpair,2); % 1: z=w1x, 2: 1/p(1-(1-p)e^(-w2x)-y)

figure()
hold on
for i=1:mt7_numpair
    yval = pairfit_info{i,2};
    xval = pairfit_info{i,3};
    zval = hecoeff_w1f(i)*xval;
    ez = (1/hecoeff_pf(i))*(1-(1-hecoeff_pf(i))*exp(-hecoeff_w2f(i)*xval)-yval);
    
    mt7_expz{i,1} = zval;
    mt7_expz{i,2} = ez;
    
    plot(exp(-zval),ez,'.','MarkerSize',12) % plot the rearranged CDF against e^z
end
plot([0,1],[0,1],'k--','LineWidth',2)
xlabel('e^{-z}')
ylabel('(1/p)(1-(1-p)e^{-(w_{2}/w_{1})z}-y)')
%ylabel('f(z,p,w_{1},w_{2})')
%ylim([-1,2])
xlim([0 1.2])
hold off





%% Check the outlier pair - is it a single exponential?
% In some datasets, mean interaction time tau goes upto a very large value
% due to artifact caused by that the pair interaction time distribution is better
% expressed by an exponential distribution than a hyperexponential
% distribution (Supplementary Info). Here we check whether the outlier pair indeed shows an
% exponential pair interaction time distribution.

% load hefit, processed_data
% here we find the bee pairs with outlier data points

smallw1_idx = find(hecoeff_w1f<0.001); % index at which tau1 > 1000
smallw2_idx = find(hecoeff_w2f<0.001); % index at which tau2 > 1000

outlierbee1_w1 = hebee1f(smallw1_idx); % bees with tau1 > 1000
outlierbee2_w1 = hebee2f(smallw1_idx); % bees with tau1 > 1000
outlierbee1_w2 = hebee1f(smallw2_idx); % bees with tau2 > 1000
outlierbee2_w2 = hebee2f(smallw2_idx); % bees with tau2 > 1000

outlierbee1 = [outlierbee1_w1; outlierbee1_w2];
outlierbee2 = [outlierbee2_w1; outlierbee2_w2];
outlierbee = [outlierbee1 outlierbee2]; % list of pairs with outlier points

beepairlist = [bee1s bee2s];
pairs = unique(beepairlist,'rows'); % unique pairs of bees
outlier_idx = zeros(size(outlierbee,1),1); % will record index of outlier bee pairs

for i=1:size(outlierbee,1)
    [~,pairidx] = ismember(outlierbee(i,:),pairs,'rows');
    outlier_idx(i) = pairidx;
end



%% fit exponential
% check whether the outliers have exponential pair interaction time
% distribution by fitting

numoutlier = size(outlierbee,1); % number of outlier pairs
r2outlier = zeros(numoutlier,1); % will record R^2 values of each fitting
toofew = 0;
expcoeff_w = zeros(numoutlier,1); % parameter for exponential, which is the rate w
expbee1 = zeros(numoutlier,1); % record whose interaction I am analyzing
expbee2 = zeros(numoutlier,1);

for i=1:numoutlier
    outlierint = fptk{outlier_idx(i),2}; % list of interaction times of the outlier pair
    [outlierecdf,xout]=ecdf(outlierint); % ecdf of outlier pair interaction times
    outlierecdf(outlierecdf==0)=[];
    xout = xout(2:end);
    
    if length(outlierecdf)>7 % minimum sample size for fitting is 3 -> 7 for better fitting
        ftexp = fittype('1-exp(-w*x)'); % CDF for exponential distribution
        expop = fitoptions(ftexp);
        expop.Lower = 0; % w>0 no negative w
        expop.StartPoint = 1;
        [expfit,expgof]=fit(xout,outlierecdf,ftexp,expop);
        
        expcoeff_w(i) = coeffvalues(expfit);
        
        beepair = fptk{outlier_idx(i),1};
        expbee1(i)=beepair(1);
        expbee2(i)=beepair(2);
    
        r2outlier(i) = expgof.rsquare;
    else
        toofew = toofew+1;
    end
end
figure()
histogram(r2outlier)
hold on
xlabel('R^{2}')
ylabel('Counts')
hold off
%% What was the hyperexponential fitting information for the outliers?
hebee = [hebee1f hebee2f];
outlier_idx_he = zeros(numoutlier,1); % index of outliers among the hyperexponential fit
outlier_g = zeros(numoutlier,1); % weight g in the hyperexponential fit
outlier_r2 = zeros(numoutlier,1); % R^2 value of the hyperexponential fit
for i=1:numoutlier
    [~,heoutidx] = ismember(outlierbee(i,:),hebee,'rows');
    outlier_idx_he(i) = heoutidx;
    outlier_g(i) = hecoeff_pf(heoutidx);
    outlier_r2(i) = r2pair_f(heoutidx);
end
%%
save('datasetName_outlierfit_date.mat','expbee1','expbee2','expcoeff_w','expfit','expgof','outlier_g','outlier_idx','outlier_idx_mt7','outlier_r2','outlierbee','r2outlier')

