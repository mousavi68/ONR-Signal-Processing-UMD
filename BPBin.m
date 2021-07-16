function [bpBin, pttBin] = BPBin(BP,PTT)
minBP = floor(min(BP));
maxBP = floor(max(BP));
bpTemp = transpose(minBP:maxBP);
pttTemp = nan(length(PTT),length(bpTemp));
for ii = 1:length(bpTemp)-1
    pttTemp(1:sum(and(BP>=bpTemp(ii),BP<bpTemp(ii+1))),ii) = PTT(and(BP>=bpTemp(ii),BP<bpTemp(ii+1)));
end
ii = length(bpTemp);
pttTemp(1:sum(BP>=bpTemp(ii)),ii) = PTT(BP>=bpTemp(ii));
pttVal = median(pttTemp,'omitnan');
bpBin = [];
pttBin = [];
for ii = 1:length(bpTemp)-1
    for jj = 1:sum(and(BP>=bpTemp(ii),BP<bpTemp(ii+1)))
        bpBin = [bpBin;bpTemp(ii)];
        pttBin = [pttBin;pttVal(ii)];
    end
end
ii = length(bpTemp);
for jj = 1:sum(BP>=bpTemp(ii))
        bpBin = [bpBin;bpTemp(ii)];
        pttBin = [pttBin;pttVal(ii)];
end
pttBin = transpose(pttBin);
end