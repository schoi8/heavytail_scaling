% get rid of last 3 digits of times, which correspond to millisec
% the camera is taking one photo per second, so resolution is 1 second.
startt2 = floor(startt/1000);
endt2 = floor(endt/1000);

% sort data according to time
roundeddata = [bee1 bee2 startt2 endt2];
rdatasort = sortrows(roundeddata,[3,4]);

bee1rs = rdatasort(:,1);
bee2rs = rdatasort(:,2);
startt2rs = rdatasort(:,3);
endt2rs = rdatasort(:,4);

% exclude data for start time > 1471928400 for 800-1, 1472187600 for 800-2
% these times are the treatment time.
%texcl = 1471928400; % 800-1
texcl = 1472187600; % 800-2
bee1e = bee1rs(startt2rs > texcl);
bee2e = bee2rs(startt2rs > texcl);
startte = startt2rs(startt2rs > texcl);
endte = endt2rs(startt2rs > texcl);
duratione = endte - startte;

% discard duration = 1. 1 second is too short to be meaningful.
bee1e2 = bee1e(duratione > 1);
bee2e2 = bee2e(duratione > 1);
startte2 = startte(duratione > 1);
endte2 = endte(duratione > 1);

% sort
exdata = [bee1e2 bee2e2 startte2 endte2];
sdata = sortrows(exdata);
bee1s = sdata(:,1);
bee2s = sdata(:,2);
startts = sdata(:,3);
endts = sdata(:,4);
durations = endts - startts;

%save('datasetName_date.mat','bee1s','bee2s','startts','endts','durations')
