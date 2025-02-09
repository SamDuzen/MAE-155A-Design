clear all; close all; clc;

%% Adjustable Parameters

Wing_Parameters = csvread("Wing_Parameters.csv");
    S = Wing_Parameters(1);
    b = Wing_Parameters(2);
    C_Root = Wing_Parameters(3);
    C_Tip = Wing_Parameters(4);
    MAC = Wing_Parameters(5);
    Y_MAC = Wing_Parameters(6);
    Sw_LE = Wing_Parameters(7);
    Sw_QC = Wing_Parameters(8);
    TaR = Wing_Parameters(9);
    ThR = Wing_Parameters(10);
    AR = Wing_Parameters(11);
    WTA = Wing_Parameters(12);
    WIA = Wing_Parameters(13);
    DA = Wing_Parameters(14);