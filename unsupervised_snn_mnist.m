% based on spnet.m code by Eugene M.Izhikevich


clearvars
close all;

% training.images(:,:,1)*255

%% clear variables and close figs
folderName = 'trial_1';
rng('default');

%% load Parameters
parameters;


%% load training set
load('mnist.mat');

%% make connection info
if isfile(strcat(folderName, '/CNX.mat'))
    load(strcat(folderName, '/CNX.mat'));
else
    CNX = connections(PARAM);
    save(strcat(folderName, '/CNX.mat'), 'CNX');

end

%% training
totalBatches = 1;

CNX.W = zeros(length(CNX.weightIdxs), totalBatches+1);
CNX.W(:, 1) = CNX.s(CNX.weightIdxs);

if ~isfolder(folderName)
    mkdir(folderName);
end

fprintf("~~~~~~~~~Starting simulation:~~~~~~~~~\n");

fileID = fopen(strcat(folderName, '/params.txt'),'w');
fprintf(fileID,'%s\n',  evalc('display(PARAM)') );
fclose(fileID);


for batch = 1:totalBatches
    
    fprintf("\tbatch: %d\n", batch);
    randSeqs = randsample(60000, 1000)';
    c = 1;
    for seqNum = randSeqs   %1:numOfSeqs %
        fprintf("\t\tseq: %d\t(%d)\n", seqNum, c);
        c = c+1;
        [CNX, cellsThatFired] = getSNN_Output(PARAM, CNX, training, seqNum);
        
        CNX.W(:, batch+1) = CNX.s(CNX.weightIdxs);
        CNX.sd = zeros(PARAM.N, PARAM.N);

    end
        
end


%% after training
save(strcat(folderName, '/PARAM.mat'), 'PARAM');
save(strcat(folderName, '/CNX_learned.mat'), 'CNX');


fprintf("finished training: %s\n", folderName);


%% determine labels in H1 layer
classCount = zeros(10, PARAM.numOfH1Cells);

for seqNum = randSeqs   %1:numOfSeqs %
    
    [~, cellsThatFired] = getSNN_Output(PARAM, CNX, training, seqNum);
    for cellNum = cellsThatFired.cellNum
        if cellNum > 28*28 && cellNum <= PARAM.numOfInputCells + PARAM.numOfH1Cells
            label = training.labels(seqNum);
            classCount(label, cellNum) = classCount(label, cellNum) + 1;
        end
    end
   
end

classification = zeros(1, PARAM.numOfH1Cells);

for i = 1: PARAM.numOfH1Cells
    [val, row] = max(classCount(:, i));
    classification(1, i) = row;
end

%% evaluating


fprintf("\tbatch: %d\n", batch);
randSeqs = randsample(60000, 60000)';
errCount = 0;

for seqNum = randSeqs(1:100)   %1:numOfSeqs %
    
    [~, cellsThatFired] = getSNN_Output(PARAM, CNX, training, seqNum);
    
    vals = zeros(1, 10);%[0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
    label = training.labels(seqNum);
    for cellNum = cellsThatFired.cellNum
        if cellNum > 28*28 && cellNum <= PARAM.numOfInputCells + PARAM.numOfH1Cells
            cellClass = classification(cellNum);
            vals(1, cellClass) = vals(1, cellClass) + 1;
        end
    end
    
    [v, col] = max(vals);
    
    if col-1 ~= label
        errCount = errCount + 1;
    end
    
end

fprintf("error count:  %d\n", errCount);



