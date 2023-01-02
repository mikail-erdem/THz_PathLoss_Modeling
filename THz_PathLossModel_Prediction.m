clc;clear;
% Specify any walue of frqeuncy: 0.7902-0.9097,0.93-0.9397 THz, 
% distance: 1:100m, height: 0:1000m
f = 0.8; % In THZ
h = 0;   % In meters
d = 6.8; % In meters
transPred = ModelTransmittance(f,h,d);

disp(['Trans-Pred Value: ',num2str(transPred)])

%% Calculate path loss
% Transmittance
[PLspread_dB_Pred , PLabs_dB_Pred] = Trans2PL(f,d,transPred);
disp(['PL-Pred Value: ',num2str(PLabs_dB_Pred)])

%% Calculate Transmittance from Model
function transPred = ModelTransmittance(f,h,d)
% Input:
% f: Frequency 0.79025-0.90974,0.93002-0.93969 THz
% h: Height 0:1000m
% d: Distance 1:100m
if ((f<0.790248000109568) || (f>0.939696000516415))
    error('The input frequency is out of modeled band {f: 0.7902-0.9097,0.93-0.9397 THz}')
end
if ((f>0.909744000434876) && (f<0.930024000490084))
    error('The input frequency is out of modeled band {f: 0.7902-0.9097,0.93-0.9397 THz}')
end
if (d<1 || d>100)
    error('The input distance is out of modeled range d: 1-100m')
end
if (h<0 || h>1000)
    error('The input height is out of modeled range d: 0-1000m')
end
% First for a given frequency "f" calculate coefficients
load Bands.mat
% BAND-1 : 0.7902-0.9097,
% BAND-2 : 0.93-0.9397
% poly8 : BAND-1
% [a1,a2,a3,a4,a5,a6,a7,a8,a9].*x.^[8,7,6,5,4,3,2,1,0]
CoeffExpFit_Band1_B1_A1=[-78091746.6392514;526711033.070705;-1553926236.73309;2619168573.80666;-2758609352.71924;1859136533.32764;-782938256.675700;188374328.578433;-19824900.2919532];
CoeffExpFit_Band1_B1_B1 = [-0.000411249080709157];
% poly4 : BAND-2
% [a1,a2,a3,a4,a5].*x.^[4,3,2,1,0]
CoeffExpFit_Band1_B2_A2 = [-74680.2770116565;280855.548765697;-396103.832267916;248296.029366205;-58368.5197872873];
CoeffExpFit_Band1_B2_B2 = [-0.000412570054095517];
% compute Coefficients
% BAND-1
[~,Idxmin] = min([min(abs(Band1-f)),min(abs(Band2-f))]);
if Idxmin == 1% f lies in BAND-1
    a = sum(CoeffExpFit_Band1_B1_A1'.*(f.^[8,7,6,5,4,3,2,1,0]));
    b = CoeffExpFit_Band1_B1_B1;
elseif Idxmin == 2% f lies in BAND-2
    a = sum(CoeffExpFit_Band1_B2_A2'.*(f.^[4,3,2,1,0]));
    b = CoeffExpFit_Band1_B2_B2;
else
    error('something wrong in the code')
end
% Now use the computed coefficent "b","a" in equation bH = a*exp(-b*x),aH=1
% For BAND-1, and BAND-2
bH = a*exp(b*h);
aH=1;
transPred = aH*exp(bH*d);
end
%% Trans-2-PathLoss
function [PLspread_dB , PLabs_dB] = Trans2PL(Freq_THz,d,Trans)
Freq_Hz = Freq_THz.*1e12;
% Trans(Trans<0)=0;
if ~isempty(Trans(Trans<0))
    error('Transmittance is less than zero at this settings')
end
SR = d;
% Absorption Coefficients K(f)
AbsCoeff = (-log(Trans))./SR;
% Total Absorption Loss [dB]
AbsorptionLoss_dB = AbsCoeff.*SR.*10.*log10(2.71828);
PLabs_dB = AbsorptionLoss_dB;
% Plspread = (SR*Freq_Hz./1e6)/10^(27.55/20);
PLspread_dB = 20.*log10(SR) + ...
    20.*log10(Freq_Hz./1e6)+20*log10(4*pi*1e6/3e8);%-GtdB-GrdB;
% Total Path Loss [dB]
end


