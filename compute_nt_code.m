%% Section nt_coding
filetext2 = fileread('nt_coding.txt');

resultLetters = regexprep(filetext2, '\s+','');

codingNtVals = nt2int(resultLetters);

% A = 1
% C = 2
% G = 3
% T = 4
codingNtVals = double(codingNtVals);

lenCode = length(codingNtVals);

% Change it to:
% A = 1
% G = 2
% C = 3
% T = 4
% because later, we have that the columns for the counts of A, G, C, T
% are in that order
for i=1:lenCode;
    if(codingNtVals(i) == 2);
        codingNtVals(i) = 3;
    elseif(codingNtVals(i) == 3);
        codingNtVals(i) = 2;
    end;
end;

nt_coding = transpose(codingNtVals);