function [nNegVals, locNegVals, negMask] = analyzeNegVals(Recon)

nNegVals = 0;
negMask = zeros(size(Recon));
for i = 1:size(Recon,1)
    for j = 1:size(Recon,2)
        if Recon(i,j)<0
            nNegVals = nNegVals + 1;
            locNegVals(1,nNegVals) = i;
            locNegVals(2,nNegVals) = j;
            negMask(i,j) = 1;
        end
    end
end