%% Main Script for force models

clear all
close all  % close all plot windows
clc  % clear command line

%% Comparison between A and B
plotPosForceModel('A', 'B', '(B-A)')
getPerturbationPosVel('A', 'B', '(B-A)')

%% Comparison between B and C
plotPosForceModel('B', 'C', '(C-B)')

%% Comparison between K and L
plotPosForceModel('K', 'L', '(L-K)')

%% Comparison between L and M
plotPosForceModel('L', 'M', '(M-L)')

%% Comparison between A and U
plotPosForceModel('A', 'U', '(U-A)')

%% Comparison between A and V
plotPosForceModel('A', 'V', '(V-A)')

%% Comparison between A and D
plotPosForceModel('A', 'D', '(D-A)')

%% Comparison between A and E
plotPosForceModel('A', 'E', '(E-A)')

%% Comparison between A and F
plotPosForceModel('A', 'F', '(F-A)')

%% Comparison between B and G
plotPosForceModel('B', 'G', '(G-B)')

%% Comparison between B and H
plotPosForceModel('B', 'H', '(H-B)')

%% Comparison between C and I
plotPosForceModel('C', 'I', '(I-C)')

%% Comparison between I and J
plotPosForceModel('I', 'J', '(J-I)')

%% Comparison between K and P
plotPosForceModel('K', 'P', '(P-K)')

%% Comparison between L and Q
plotPosForceModel('L', 'Q', '(Q-L)')

%% Comparison between L and R
plotPosForceModel('L', 'R', '(R-L)')

%% Comparison between M and S
plotPosForceModel('M', 'S', '(S-M)')

%% Comparison between S and T
plotPosForceModel('S', 'T', '(T-S)')

%% Comparison between D and N
plotAccForceModel('D', 'N', '(N-D)')
getPerturbationAcc('D', 'N', '(N-D)')

%% Comparison between K and N
plotAccForceModel('N', '', 'N (Atmosferic Drag)')

%% Comparison between D and O
plotAccForceModel('E', 'O', '(O-E)')

%% Comparison between K and O
plotAccForceModel('O', '', 'O (Solar Radiation Pressure)')

%% Air density N
plotDensityForceModel('N', 'N')