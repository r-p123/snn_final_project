% based on spnet.m code by Eugene M.Izhikevich

function CNX = connections(PARAM)

CNX.delays = cell(PARAM.N, PARAM.D);
CNX.post = zeros(PARAM.N, PARAM.N);


CNX.s = -100*ones(PARAM.N, PARAM.N);

%% connections from input neurons to hidden layer


% CNX.delay_vals = -1*ones(PARAM.N, PARAM.N);


% connect input to H1 layer
for i = PARAM.numOfInputCells + 1 : PARAM.numOfH1Cells + PARAM.numOfInputCells
    
    for c = 1:PARAM.numOfInputCells
        CNX.post(c, i) = i;
        CNX.s(c, i) = PARAM.initialWeight;
        
        CNX.delays{c, 1}(end+1) = i;
        CNX.delay_vals(c, i) = 0;
    end
    
end

clear x c i


% connect H1 to H2 layer
c = PARAM.numOfInputCells + 1;
for i = PARAM.numOfInputCells + PARAM.numOfH1Cells + 1 : PARAM.numOfH1Cells + PARAM.numOfH2Cells + PARAM.numOfInputCells
    
    CNX.post(c, i) = i;
    CNX.s(c, i) = PARAM.initialWeight;
    
    CNX.delays{c, 1}(end+1) = i;
    CNX.delay_vals(c, i) = 0;
        
    c = c+1;
    
end

% connect inhibitory layer to H1 layer
for c = PARAM.numOfInputCells + PARAM.numOfH1Cells + 1 : PARAM.numOfH1Cells + PARAM.numOfH2Cells + PARAM.numOfInputCells
    
    for i = PARAM.numOfInputCells + 1 : PARAM.numOfInputCells  + PARAM.numOfH1Cells
        CNX.post(c, i) = i;
        CNX.s(c, i) = -1*PARAM.initialWeight;
        
        CNX.delays{c, 1}(end+1) = i;
        CNX.delay_vals(c, i) = 0;
    end
    
end

%% weights

CNX.sd=zeros(PARAM.N, PARAM.N);

CNX.weightIdxs = [];

for i = 1: PARAM.N * PARAM.N
    if CNX.s(i) >=0
        CNX.weightIdxs(end+1) = i;
    end
end

clear i

%% information for STDP weight changes

CNX.pre = cell(PARAM.N,1);
CNX.aux = cell(PARAM.N,1);

for i=1 : PARAM.N
    for j=1 : PARAM.D
        for k=1:length(CNX.delays{i,j})
            CNX.pre{ CNX.delays{i, j}(k) }(end+1) = PARAM.N * (CNX.delays{i, j}(k)-1)+i;
            CNX.aux{ CNX.delays{i, j}(k) }(end+1) = PARAM.N * (PARAM.D-1-j)+i; % takes into account delay
        end
    end
end

clear i j k


end