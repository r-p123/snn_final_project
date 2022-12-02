%parameters
% based on spnet.m code by Eugene M.Izhikevich

PARAM.noiseFlag = 0;

PARAM.initialWeight = 0.50; %0.76 IT STARTS SPIKING



PARAM.resetVoltage = -2;
PARAM.tau_m = 10;
PARAM.tau_s = PARAM.tau_m/4;
PARAM.K1 = 2;
PARAM.K2 = 4;

PARAM.thresh = 2;


PARAM.Tmax = 350;  %100;   %
PARAM.dt = 1;
PARAM.timeVals = 0:PARAM.dt:PARAM.Tmax;
PARAM.length_t = length(PARAM.timeVals);


PARAM.numOfH1Cells = 400;
PARAM.numOfH2Cells = 400;
PARAM.numOfInputCells = 784;
PARAM.N = PARAM.numOfH1Cells + PARAM.numOfH2Cells + PARAM.numOfInputCells;

PARAM.D = 1;
PARAM.maxW = 1;
PARAM.decayRate = 0;


PARAM.a_plus = 0.03; %0.03;
PARAM.a_minus = 1.2*PARAM.a_plus;  %0.85 * param.a_plus;

PARAM.LTP_Values = PARAM.a_plus * exp(-1*PARAM.timeVals(1:100)/10);  %exp(-1*param.timeVals(1:50)/5);
PARAM.LTD_Values = PARAM.a_minus * exp(-1*PARAM.timeVals(1:100)/30);  %exp(-1*param.timeVals(1:50)/10);


PARAM.STDP = 1;
% PARAM.decayRate = 0.0001;


%% calculate v0
tVals = 0:1:100;
kernelPSP = zeros(1, length(tVals));

for i = 1 : length(tVals)   
    epsilon = ( exp(-(i)/PARAM.tau_m) -  exp(-(i)/PARAM.tau_s) );
    kernelPSP(i) = epsilon;  
end

PARAM.v0 = 1/max(kernelPSP);

clear tVals kernelPSP epsilon i
