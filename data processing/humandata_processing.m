%% collect the duration instead of the exact time the interaction happened
% for human data that provide the time of interaction detection rather than
% duration

% variables
% bee1: id for first agent
% bee2: id for second agent
% time: the time at which interaction was detected

rawdata = [bee1 bee2 time];
beedatasort = sortrows(rawdata); 
times = beedatasort(:,3);

n = size(bee1,1);
pair = beedatasort(1,1:2); % first pair in the list
beedatasort2 = zeros(n,3); % this will be the new matrix of the form [bee1 bee2 duration]
pairind = 1; % index for a pair interaction
dt = 20; % duration of interaction. In SocioPatterns recording, the time step is 20 seconds.
for i=2:n
    if beedatasort(i,1:2)==pair % if the pair is the same
        if times(i)==times(i-1)+20 % still the same interaction
            dt = dt+20;
        else % different interaction event
            beedatasort2(pairind,3)=dt;
            dt = 20;
            beedatasort2(pairind,1:2) = pair;
            pairind = pairind+1;
        end
    else % different pair
        beedatasort2(pairind,3)=dt;
        dt = 20;
        beedatasort2(pairind,1:2) = pair;
        pair = beedatasort(i,1:2);
        pairind = pairind+1;
    end
end
dummyvec = beedatasort2(:,1);
beedatasort2(dummyvec==0,:)=[]; % get rid of rows with 0 entries

bee1s = beedatasort2(:,1); % new sorted column representing agent 1
bee2s = beedatasort2(:,2); % new sorted column representing agent 2
durations = beedatasort2(:,3); % new sorted column representing duration

