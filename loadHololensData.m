function [vicon, pvhololens, alldata, allhololens] = loadHololensData(folder, viconPath)
% loads all the HoloLens data from the folder and vicon data from viconPath
    vicon = readtable(viconPath);
    pvhololens = readtable([folder, 'pv.csv']);
    vlcll = readtable([folder, 'vlc_ll.csv']);
    vlclf = readtable([folder, 'vlc_lf.csv']);
    vlcrf = readtable([folder, 'vlc_rf.csv']);
    vlcrr = readtable([folder, 'vlc_rr.csv']);
    longthrowdepth = readtable([folder, 'long_throw_depth.csv']);

    % all cameras from HoloLens in one array
    allhololens = [pvhololens; vlcll; vlclf; vlcrf; vlcrr; longthrowdepth];
    allhololens = sortrows(allhololens, 'Timestamp');

    % cell array of each HoloLens camera
    alldata = {pvhololens, vlcll, vlclf, vlcrf, vlcrr, longthrowdepth};

end
