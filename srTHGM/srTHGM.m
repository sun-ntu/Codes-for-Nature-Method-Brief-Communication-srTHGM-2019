clear all;
close all;
clc;

Foldername = pwd;
demo = input('Demo data? (1) mice (2) human (3) all :','s');
if(demo=='1')
    subFolder = 'mice';
    srTHGM_calibration;
    srTHGM_imageprocessing;
elseif(demo=='2')
    subFolder = 'human';
    srTHGM_calibration;
    srTHGM_imageprocessing;
elseif(demo=='3')
    subFolder = 'mice';
    srTHGM_calibration;
    srTHGM_imageprocessing;
    subFolder = 'human';
    srTHGM_calibration;
    srTHGM_imageprocessing;
else
    disp('Wrong selection! Please run program again.');
end