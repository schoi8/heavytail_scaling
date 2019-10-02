%% For unsorted bee data from the 2013 experiment
% process the raw data to be analyzed

% variables
% bee1 = id1, bee2 = id2, startt = begin, endt = end

endt2 = floor(endt/1000);
startt2 = floor(startt/1000);
% divide by 1000 because the resolution of the camera is 1 second so the
% millisecond times in the unix epoch timestamp are not meaningful. 
rawdata = [bee1 bee2 startt2 endt2];
bdatasort = sortrows(rawdata);
bee1st = bdatasort(:,1);
bee2st = bdatasort(:,2);
starttst = bdatasort(:,3);
endtst = bdatasort(:,4);
durationst = endtst - starttst;

% discard duration = 1. 1 second is too short to be meaningful.
bee1s = bee1st(durationst > 1);
bee2s = bee2st(durationst > 1);
startts = starttst(durationst > 1);
endts = endtst(durationst > 1);
durations = durationst(durationst > 1);

%save('datasetName_sorted_date.mat','bee1s','bee2s','startts','endts','durations')
