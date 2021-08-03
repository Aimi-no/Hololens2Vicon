function res = optimizeAlignmentTuning(Binit, Bmax, pcVicon, ViconRot, allPointcloudData)
% init
res.minB = 0;
res.params = [];
res.vals = [];
res.minval = realmax;
res.initvals = [];
initpar = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0];

% iterate over all possible B
for B = Binit:Bmax 
    % get Vicon markers and rotations dependent on HoloLens timestamp steps
    vic = cell(6,1);
    Rvic = cell(6,1);
    vicsize = size(pcVicon.Location, 1);
    for k = 1:6
        i = allPointcloudData.cs{k} + B; 
        i(i > vicsize) = [];
        vic{k} = pcVicon.Location(i, :);
        Rvic{k} = ViconRot(i, :);
    end
    
    % init optimization
    % parameters: t, rho, S, d, C, R, D
    fun = @(par)closerTransformation([[par(1) par(2) par(3)]', [par(4) par(5) par(6)]', ...
    [par(7) par(8) par(9)]', [par(10) par(11) par(12)]', [par(13) par(14) par(15)]', ...
    [par(16) par(17) par(18)]'], par(19), euler2mat([par(20) par(21) par(22)]), [par(23) par(24) par(25)]', ...
    vic, Rvic, allPointcloudData.hol);

    options = optimset( 'MaxFunEvals', 1000 * 25);
    % optimize
    [bestpar, bestval] = fminsearch(fun, initpar, options);
    

    % save all results
    res.initvals = [res.initvals closerTransformation([[initpar(1) initpar(2) initpar(3)]', [initpar(4) initpar(5) initpar(6)]', ...
        [initpar(7) initpar(8) initpar(9)]', [initpar(10) initpar(11) initpar(12)]', [initpar(13) initpar(14) initpar(15)]', ...
        [initpar(16) initpar(17) initpar(18)]'], initpar(19), euler2mat([initpar(20) initpar(21) initpar(22)]), [initpar(23) initpar(24) initpar(25)]', vic, Rvic, allPointcloudData.hol)];
    
    res.params = [res.params bestpar];
    res.vals = [res.vals bestval];
    
    % save best results
    if bestval < res.minval
      minpar = bestpar;
      res.minval = bestval; 
      res.minB = B;
    end
    fprintf(['For B = ', num2str(B), ' the optimized value is ',  num2str(bestval),' parameters are:', '\n' ]);
    bestpar
    
end

fprintf(['Best B is ', num2str(res.minB), ' the optimized value is ',  num2str(res.minval),' parameters are:', '\n' ]);
minpar

% parse optimized parameters
res.St = euler2mat([minpar(20) minpar(21) minpar(22)])';
res.rho = minpar(19);
res.d = [minpar(23) minpar(24) minpar(25)]';
res.t = [[minpar(1) minpar(2) minpar(3)]' ...
[minpar(4) minpar(5) minpar(6)]' ...
[minpar(7) minpar(8) minpar(9)]' ...
[minpar(10) minpar(11) minpar(12)]' ...
[minpar(13) minpar(14) minpar(15)]' ...
[minpar(16) minpar(17) minpar(18)]'];

end

function res = closerTransformation(t, rho, S, d, C, Rv, D)
% calculates parameters for alignment tuning
%     C - known, Vicon marker
%     R(e) - known, Vicon rotation
%     t - unknown, translation between Vicon marker and HoloLems camera 
%     rho - unknown, scale of the Viconvectors onto HoloLens vectors
%     S(e) - unknown, camera rotation in HoloLens coordinate system
%     D - known, HoloLens camera centers
%     d - unknown, origin of HoloLens coordinate system in Vicon coordinate system

    res = cell(6,1);

    % cycle over all cameras
    for k = 1:6

        
        num = size(Rv{k},1);
        res{k} = zeros(num, 1);

        for i = 1:num
            R = euler2mat(Rv{k}(i,:));
            % objective function
            f = C{k}(i,:)' + R * t(:,k) - (1/rho) * S' * D{k}(i,:)' - d;
            % errors
            res{k}(i) = norm(f)^2;
        end
    end
    res = cell2mat(res);
    
    % filtering out outliers caused by Vicon errors
    idx = kmeans(res, 2, 'Start', [min(res); max(res)]);
    idxmin = idx(res == min(res));
    
    % final error
    res = sum(res(idx == idxmin(1))) + 10 * (norm(t) + norm(d));
    
end