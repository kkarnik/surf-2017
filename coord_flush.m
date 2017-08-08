exonCoords = dlmread('exonCoordData.txt');

exonCoordsWithPosns = [exonCoords zeros(size(exonCoords, 1),1)];

% First find out which direction the gene is read in
% If 1, then in TSC1 direction. If 0, then in TS2 direction.
TSC1direction = (exonCoords(1, 2) > exonCoords(2, 2));

numRows = size(exonCoords, 1);

codonIndex = 0;

if(TSC1direction == 1);
    for i=numRows-1:-1:2;
        codonIndex = codonIndex + 1;
        exonCoordsWithPosns(i, 5) = codonIndex;
        codonIndex = codonIndex + exonCoords(i, 4) - exonCoords(i, 3);
    end;
else
    for i=6:1:numRows-1;
        codonIndex = codonIndex + 1;
        exonCoordsWithPosns(i, 5) = codonIndex;
        codonIndex = codonIndex + exonCoords(i, 4) - exonCoords(i, 3);
    end;
end;