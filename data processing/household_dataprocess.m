% m1 and m2 are id's for agents in household data.
% household m1 and m2 pair not always m1<m2 -> change the order
m1order = m1;
m2order = m2;
for i=1:length(m1)
    if m1(i)>m2(i)
        m1order(i) = m2(i);
        m2order(i) = m1(i);
    end
end

%% assign new unique id's for household data so that there is no overlap of id's across households
% it turns out that same values of m1 and m2 were across different
% households so I assigned new unique id's not overlapping across
% households

% separate data according to household: L,F,E,B,H
house1 = char(h1);
house2 = char(h2);

h1L_idx = find(house1=='L');
h1F_idx = find(house1=='F');
h1E_idx = find(house1=='E');
h1B_idx = find(house1=='B');
h1H_idx = find(house1=='H');

mL1 = m1order(h1L_idx);
mF1 = m1order(h1F_idx);
mE1 = m1order(h1E_idx);
mB1 = m1order(h1B_idx);
mH1 = m1order(h1H_idx);
mL2 = m2order(h1L_idx);
mF2 = m2order(h1F_idx);
mE2 = m2order(h1E_idx);
mB2 = m2order(h1B_idx);
mH2 = m2order(h1H_idx);
durL = duration(h1L_idx);
durF = duration(h1F_idx);
durE = duration(h1E_idx);
durB = duration(h1B_idx);
durH = duration(h1H_idx);

rawL = [mL1 mL2 durL];
rawF = [mF1 mF2 durF];
rawE = [mE1 mE2 durE];
rawB = [mB1 mB2 durB];
rawH = [mH1 mH2 durH];

sortL = sortrows(rawL);
sortF = sortrows(rawF);
sortE = sortrows(rawE);
sortB = sortrows(rawB);
sortH = sortrows(rawH);

%% assign unique id for individuals not overlapping across household
msL1 = sortL(:,1);
msL2 = sortL(:,2);
uniqL = unique([msL1;msL2]);
msF1 = sortF(:,1);
msF2 = sortF(:,2);
uniqF = unique([msF1;msF2]);
msE1 = sortE(:,1);
msE2 = sortE(:,2);
uniqE = unique([msE1;msE2]);
msB1 = sortB(:,1);
msB2 = sortB(:,2);
uniqB = unique([msB1;msB2]);
msH1 = sortH(:,1);
msH2 = sortH(:,2);
uniqH = unique([msH1;msH2]);

lengL = length(msL1);
lengF = length(msF1);
lengE = length(msE1);
lengB = length(msB1);
lengH = length(msH1);

newL1 = zeros(lengL,1);
newL2 = zeros(lengL,1);
newF1 = zeros(lengF,1);
newF2 = zeros(lengF,1);
newE1 = zeros(lengE,1);
newE2 = zeros(lengE,1);
newB1 = zeros(lengB,1);
newB2 = zeros(lengB,1);
newH1 = zeros(lengH,1);
newH2 = zeros(lengH,1);
%%
numL = length(uniqL);
for i=1:numL
    findL1 = find(msL1==uniqL(i));
    findL2 = find(msL2==uniqL(i));
    if sum(findL1)>0
        newL1(findL1)=i;
    end
    if sum(findL2)>0
        newL2(findL2)=i;
    end
end

numF = length(uniqF);
for i=1:numF
    findF1 = find(msF1==uniqF(i));
    findF2 = find(msF2==uniqF(i));
    if sum(findF1)>0
        newF1(findF1)=i+numL;
    end
    if sum(findF2)>0
        newF2(findF2)=i+numL;
    end
end

numE = length(uniqE);
for i=1:numE
    findE1 = find(msE1==uniqE(i));
    findE2 = find(msE2==uniqE(i));
    if sum(findE1)>0
        newE1(findE1) = i+numL+numF;
    end
    if sum(findE2)>0
        newE2(findE2) = i+numL+numF;
    end
end

numB = length(uniqB);
for i=1:numB
    findB1 = find(msB1==uniqB(i));
    findB2 = find(msB2==uniqB(i));
    if sum(findB1)>0
        newB1(findB1) = i+numL+numF+numE;
    end
    if sum(findB2)>0
        newB2(findB2) = i+numL+numF+numE;
    end
end

numH = length(uniqH);
for i=1:numH
    findH1 = find(msH1==uniqH(i));
    findH2 = find(msH2==uniqH(i));
    if sum(findH1)>0
        newH1(findH1) = i+numL+numF+numE+numB;
    end
    if sum(findH2)>0
        newH2(findH2) = i+numL+numF+numE+numB;
    end
end

bee1 = [newL1;newF1;newE1;newB1;newH1];
bee2 = [newL2;newF2;newE2;newB2;newH2];
duration = [sortL(:,3);sortF(:,3);sortE(:,3);sortB(:,3);sortH(:,3)];
