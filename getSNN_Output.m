% based on spnet.m code by Eugene M.Izhikevich


function [CNX, cellsThatFired] = getSNN_Output(PARAM, CNX, training, seqNum)




PSP = zeros(PARAM.N, PARAM.Tmax+1+PARAM.D);

LTPot = zeros(PARAM.N, PARAM.Tmax+1+PARAM.D + 50);
LTDep = zeros(PARAM.N, PARAM.Tmax+1+PARAM.D + 50);
v = zeros(PARAM.N,1);       % initial values
v_info = cell(PARAM.N, 2);  % need presyn cell number and spike time to calculate EPSP
lastPostSynSpikeTime = (PARAM.Tmax+100)*ones(PARAM.N, 1);
refractoryUntil = -1*ones(PARAM.N, 1);

firings=[-PARAM.D 0];                         % spike timings

idx1 = 1;


A = training.images(:,:,seqNum)*255/4;
A = reshape(A,1,28*28);

inputSpikeTimes = [];
inputCellNum = [];
for i = 1 : 28*28
    fr = A(i); %fr -> Hz
    if fr>0
        dt = 1/1000; % s
        nBins = PARAM.Tmax; % 10 ms spike train
        spikeTrain = rand(1, nBins) < fr*dt;
        [~, idxs] = find(spikeTrain);
        inputSpikeTimes = [idxs, inputSpikeTimes];
        inputCellNum = [i*ones(1, length(idxs)) , inputCellNum];
    end
end
A = [inputSpikeTimes; inputCellNum ]';
B = sortrows(A, 1);
inputSpikeTimes = B(:, 1);
inputCellNum = B(:, 2);
clear A B

for t = 1 : 1 : PARAM.Tmax
    
    
    %incorporate input spike train
    while idx1 <= length(inputSpikeTimes) && inputSpikeTimes(idx1) == t
        cellNum = inputCellNum(idx1);
        v(cellNum) = PARAM.thresh;
        idx1 = idx1+1;
    end
 
    
    fired = find(v>=PARAM.thresh & v~=PARAM.thresh*PARAM.K1);                % indices of fired neurons
    
    for k = 1 : length(fired)
        if fired(k) > 28*28
            
            if PARAM.STDP
                CNX.sd(CNX.pre{fired(k)}) = CNX.sd(CNX.pre{fired(k)}) + LTPot(PARAM.N*t + CNX.aux{fired(k)});
            end
            
           
            
            v_info{fired(k), 1} = [];
            v_info{fired(k), 2} = [];
           
        end
    end
    
    v(fired) = PARAM.resetVoltage;
    lastPostSynSpikeTime(fired) = t;
    refractoryUntil(fired) = t+5;
    
    if ~isempty(fired)
        for k = 1 : length(fired)
            LTPot(fired(k), t+PARAM.D : t+PARAM.D+49 ) = PARAM.LTP_Values( 1 : 50 ); % irrelevant for cells 1-6 b/c input cells
            LTDep(fired(k), t+PARAM.D : t+PARAM.D+49 ) = PARAM.LTD_Values( 1 : 50 );
        end
    end
    
    
    firings=[  firings ; t*ones(length(fired),1) , fired  ];
    k=size(firings,1);
    
    while firings(k,1) > t-PARAM.D
        
        % firings(k,1) = time
        % firings(k,2) = cell num
        cellsDownTheLine = CNX.delays{firings(k,2), t-firings(k,1)+1 };
        ind = CNX.post(firings(k,2),cellsDownTheLine);
        
        for c = cellsDownTheLine
            
            if t > refractoryUntil(c)
                v_info{c, 1}(end+1) = t;
                v_info{c, 2}(end+1) = firings(k,2);
            end
            
        end
        
        if PARAM.STDP
            CNX.sd(firings(k,2), cellsDownTheLine) = CNX.sd(firings(k,2), cellsDownTheLine) - LTDep(ind, t+PARAM.D)';
        end
        
        
        k=k-1;
        
    end
    
    %calculate epsp for all hidden cells
    for i = 1 : PARAM.N
        
        len = length(v_info{i, 1});
        
        ti = lastPostSynSpikeTime(i);
        
        if t-ti >= 0 && i > 6
            eta = PARAM.thresh * (PARAM.K1 * exp(-(t - ti)/PARAM.tau_m) - PARAM.K2*(exp(-(t - ti)/PARAM.tau_m) - exp(-(t - ti)/PARAM.tau_s)) );
        else
            eta = 0;
        end
        
        epsTot = 0;
        for k = 1: len
            tj = v_info{i, 1}(k);
            wi = CNX.s(v_info{i, 2}(k), i);
            epsilon = wi*PARAM.v0*( exp(-(t - tj)/PARAM.tau_m) -  exp(-(t-tj)/PARAM.tau_s) );
            epsTot = epsTot + epsilon;
        end
        
        
        v(i) = eta + epsTot;
        PSP(i, t) = v(i);
    end
    clear i;
end

cellsThatFired.cellNum = firings(2:end, 2)';
cellsThatFired.spikeTime = firings(2:end, 1)';


if PARAM.STDP
    for widx = CNX.weightIdxs
        CNX.s(widx) = max(0, min(PARAM.maxW,    CNX.s(widx) + CNX.sd(widx) ) );
    end
end


end


