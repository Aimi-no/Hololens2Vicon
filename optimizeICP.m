function res = optimizeICP(pcVicon, pcHoloLens)
% align Vicon and HoloLens poses using random rotation and ICP over 100
% iterations

    res.besterr = realmax;
    res.bestax = 0;
    res.bestay = 0;
    res.bestaz = 0;
    res.bestrmse = 0;
    res.R = diag([1 1 1]);

    % build a KD-tree searcher of Vicon to measure alignment error
    MdlICP = KDTreeSearcher(pcVicon.Location);

    for i = 1:100
        % random rotation in euler angles
        ay = interp1([0,1],[0,360],rand(1));
        az = interp1([0,1],[0,360],rand(1));
        ax = interp1([0,1],[0,360],rand(1));


        Rx = [1 0 0 0; 0 cos(ax) -sin(ax) 0; 0 sin(ax) cos(ax) 0; 0 0 0 1];
        Ry = [cos(ay) 0 sin(ay) 0; 0 1 0 0; -sin(ay) 0 cos(ay) 0; 0 0 0 1];
        Rz = [cos(az) -sin(az) 0 0; sin(az) cos(az) 0 0; 0 0 1 0; 0 0 0 1];

        tmpR = Rx * Ry * Rz;

        % rotate HoloLens cameras
        tform_rotate = affine3d(tmpR);
        ptrotHolol = pctransform(pcHoloLens,tform_rotate);
        
        % run ICP
        [tform,tmphololensReg, rmse] = pcregrigid(ptrotHolol, pcVicon, 'MaxIterations', 100000, 'Tolerance', [0.001, 0.0009]);

        % NN search of HoloLens cams in Vicon to estimate ICP error
        [~, D] = knnsearch(MdlICP, tmphololensReg.Location);
        err = sum(D);

        % save the best result
        if err < res.besterr
            res.besterr = err;
            res.besttform = tform;
            res.hololensReg = tmphololensReg;
            res.bestax = ax;
            res.bestay = ay;
            res.bestaz = az;
            res.bestrmse = rmse;
            res.besttform_rotate = tform_rotate;
            res.R = tmpR;
            res.bestD = D;
        end
    end
end