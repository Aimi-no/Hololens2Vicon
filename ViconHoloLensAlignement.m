%% ---------- Init and loading of Vicon and HoloLens data ----------
croppingStart = 0;
croppingEnd = 100;
Binit = 0;
Bmax = 100;
numOfCameras = 6;
folder = './Vicon_session_2020_12_02/HoloLensRecording__2020_12_02__12_54_39/';
viconPath = './Vicon_session_2020_12_02/hololens_seq03.txt';

% load all data
[vicon, pvhololens, alldata, allhololens] = loadHololensData(folder, viconPath);

% filter out invalid poses in Vicon
indexes = vicon.Var4(:) ~= 1;
pcVicon = pointCloud([vicon.Var5(indexes), vicon.Var6(indexes), vicon.Var7(indexes)]);
ViconRot = [vicon.Var5(indexes), vicon.Var6(indexes), vicon.Var7(indexes)];
pcHoloLens = pointCloud([allhololens.Position_X(croppingStart+1:end-croppingEnd), allhololens.Position_Y(croppingStart+1:end-croppingEnd), allhololens.Position_Z(croppingStart+1:end-croppingEnd)]);

% show Vicon and HoloLens poses
figure();
pcshowpair(pcHoloLens, pcVicon, 'MarkerSize', 50);

axis equal;
xlabel("x");
ylabel("y");
zlabel("z");
title('Hololens and Vicom tracking before registration');
legend('\color{white} Hololens', '\color{white} Vicon');

%% ---------- ICP alignment ---------- 

% run ICP optimization
resICP = optimizeICP(pcVicon, pcHoloLens);

% show alignment after ICP
figure();
pcshowpair(resICP.hololensReg, pcVicon, 'MarkerSize', 50);
axis equal;
xlabel("x");
ylabel("y");
zlabel("z");
title(['Registered Hololens tracking and Vicom tracking, ax = ', num2str(resICP.bestax), ', ay = ', num2str(resICP.bestay), ', az = ', num2str(resICP.bestaz)]);
legend('\color{white} Hololens', '\color{white} Vicon');

%% ---------- Objective function optimization ---------- 

%numHolol = size(pvhololens,1);
%numVicon = size(vicon,1);

% create pointclouds and related information
allPointcloudData = createPointclouds(alldata, numOfCameras, resICP.besttform_rotate, resICP.besttform);

% optimize parameters of the objective function
optimizedPars = optimizeAlignmentTuning(Binit, Bmax, pcVicon, ViconRot, allPointcloudData);

%% ---------- Show graphs and results of objective function optimization ---------- 

% plot function values(errors)
figure();
plot(Binit:Bmax , optimizedPars.vals, 'b', Binit:Bmax , optimizedPars.initvals, 'r');
legend('Optimized function value', 'Initial function value');
xlabel('\beta');
ylabel('Error over all cameras [m]')


errs = [];
hol2viccorr = cell(6,1);

% plot tuned alignment over all cameras
for m = 1:numOfCameras
    
    % init transformed HoloLens camera to Vicon
    hol2vicon = zeros(size(allPointcloudData.hol{m},1), 3);
    
    % timestamps in Vicon corresponding to Hololens timestamps
    j = allPointcloudData.cs{m} + optimizedPars.minB; 
    j(j > pcVicon.Count) = [];
    
    % init Vicon and shifted Vicon
    vic = pcVicon.Location(j, :);
    victrans = zeros(size(allPointcloudData.hol{m},1), 3);

    Rvic = ViconRot(j, :);

    % transform HoloLens camera to Vicon
    hol2vicon = (1/optimizedPars.rho * optimizedPars.St * allPointcloudData.hol{m}' - optimizedPars.d)';
    
    % shift Vicon according to t
    for i = 1:size(allPointcloudData.hol{m},1)
        victrans(i, :) = (vic(i, :)' + euler2mat(Rvic(i, :)) * optimizedPars.t(:,m))';
    end

    pcHol = pointCloud(hol2vicon);

    % plot Vicon and HoloLens alignment with errors
    figure();
    pcshowpair(pcHol, pcVicon, 'MarkerSize', 50);
    axis equal;
    xlabel("x");
    ylabel("y");
    zlabel("z");
    title(['Registered Hololens tracking and Vicom tracking after tuning for ', num2str(m), 'th camera']);
    hold on;

    pcshow(pointCloud(allPointcloudData.hol{m}, 'Color', repmat([1 1 0], length(allPointcloudData.hol{m}),1)), 'MarkerSize', 50);
    pcshow(pointCloud(victrans, 'Color', repmat([0 1 1], length(victrans), 1)),'MarkerSize', 50);
    hold on;

    % plot errors
    for i = 1:pcHol.Count
        if(norm(victrans(i,:) - hol2vicon(i,:)) < 0.2) % filtering out outliers
            errs = [errs norm(victrans(i,:) - hol2vicon(i,:))];
            plot3([victrans(i, 1), hol2vicon(i,1)], [victrans(i, 2), hol2vicon(i,2)], [victrans(i, 3), hol2vicon(i,3)], 'r');
        end
    end

    legend('\color{white} Transformed HoloLens','\color{white} Vicon', '\color{white} Hololens before tuning', '\color{white} Shifted Vicon', '\color{white} Error between HoloLens');
    title(['Registered Hololens tracking and Vicom tracking after tuning for ', num2str(m), 'th camera']);

    % calculate errors
    dist = vecnorm((victrans - hol2vicon)');
    % filter out outliers
    distgroups = kmeans(dist', 2);
    idxmin = distgroups(dist == min(dist));
    
    % print mean error
    meanerr = mean(dist(distgroups == idxmin));
    fprintf(['For ', num2str(m), 'th camera the mean distance is: ',  num2str(meanerr), '\n' ]);
    hol2viccorr{m} = victrans - hol2vicon;

end

% plot all errors in alignment over all Hololens cameras
figure();
plot(sort(errs));
xlabel("err id");
ylabel("error[m]");
title('All errors between aligned HoloLens and Vicon');

%% ---------- Transform depth data ---------- 

% load all depth files
files = dir([folder, 'long_throw_depth/']);
depthposes = readtable([folder, 'long_throw_depth.csv']);
depthposes(1,:) = [];

% load Matterport pointcloud
pcMatterport = pcread('D:/Documents/Work/B-315/matterport2vicon.ply');
pointclouds_depth = {};

ind = 1;
cols = getColors(length(files));

% get Vicon poses corresponding to the best B
j = allPointcloudData.cs{6} + optimizedPars.minB; 
Rvic = ViconRot(j, :);
vic = pcVicon.Location(j, :);

% get rotations used in ICP registration together
Rii = resICP.R(1:3, 1:3) * resICP.besttform.Rotation;

% get all translations together
T = optimizedPars.St * resICP.besttform.Translation' - optimizedPars.d;

% all rotation between HoloLens and Vicon coordinate systems
allqs = [];
%figure();

% over all files in folder
for i = 1:length(files)
    % selecting correct files
    if ~startsWith(files(i).name, 'world_') && endsWith(files(i).name, '.ply') && ~endsWith(files(i).name,'_transf.ply')
        % parsing name
        name = files(i).name;
        parsed = split(name, '.');
        poseID = ['long_throw_depth\' , parsed{1}, '.pgm'];
        % getting information about the current file from the csv depth file
        row = depthposes(strcmp(depthposes.ImageFileName, poseID), :);
        
        % if file is in the table
        if ~isempty(row)
            % get the serial number of the pose (so we can use coresponding Vicon)
            [~, rowindex] = ismember(poseID, depthposes.ImageFileName);
            
            % get HoloLens pose matrices
            FrameToOrigin = [row.FrameToOrigin_m11 row.FrameToOrigin_m12 row.FrameToOrigin_m13 row.FrameToOrigin_m14;
                    row.FrameToOrigin_m21 row.FrameToOrigin_m22 row.FrameToOrigin_m23 row.FrameToOrigin_m24;
                    row.FrameToOrigin_m31 row.FrameToOrigin_m32 row.FrameToOrigin_m33 row.FrameToOrigin_m34;
                    row.FrameToOrigin_m41 row.FrameToOrigin_m42 row.FrameToOrigin_m43 row.FrameToOrigin_m44];

            CameraViewTransform = [row.CameraViewTransform_m11 row.CameraViewTransform_m12 row.CameraViewTransform_m13 row.CameraViewTransform_m14;
                          row.CameraViewTransform_m21 row.CameraViewTransform_m22 row.CameraViewTransform_m23 row.CameraViewTransform_m24;
                          row.CameraViewTransform_m31 row.CameraViewTransform_m32 row.CameraViewTransform_m33 row.CameraViewTransform_m34;
                          row.CameraViewTransform_m41 row.CameraViewTransform_m42 row.CameraViewTransform_m43 row.CameraViewTransform_m44];


            C2D = inv(CameraViewTransform)';
            D2C = CameraViewTransform';
            D2O = FrameToOrigin'; 
            O2D = inv(D2O);
            Rhi = D2O * C2D;

            oldRhi = Rhi(1:3, 1:3);
            mChw = oldRhi' * Rhi(1:3,4);

            % load HoloLens depth pointcloud
            pcdepth = pcread([folder, 'long_throw_depth/', name]);
            depthdata = pcdepth.Location / -1000;
            pcdepth = pointCloud(depthdata);
            corrections = hol2viccorr{6};

            % filtering out outliers
            if(norm(corrections(rowindex,:)) < 0.2) 
                % calculating rotation between HoloLens and Vicon
                % coordinate systems
                
                Rv = euler2mat(Rvic(rowindex,:));
                Rhv = optimizedPars.St * Rii' * oldRhi;
                Rdelta = Rhv * Rv';
                qdelta = r2q(Rdelta); 
                allqs = [allqs qdelta];
                allRs{ind} =  Rdelta;

                Rhv_star = Rdelta * Rv';

                Rstar_test = Rv' * optimizedPars.St * Rii' * oldRhi';
                Rstarq{ind} = r2q(Rstar_test);
                % we want Rhv = Rdelta * Rv' to hold;
                
                % transforming HoloLens depth pointcloud to Vicon 
                scale = 1/optimizedPars.rho;
                rotation = optimizedPars.St * Rii';
                translation = T + Rv * optimizedPars.t(:,6) + corrections(rowindex,:)';
                depthdata = transformDepthDataToVicon(depthdata, D2O * C2D, scale, rotation, translation);
                
                % save transformed pointcloud
                cmatrix = ones(pcdepth.Count, 1) * cols{i};   
                pcDepthTransformed = pointCloud(depthdata(1:1000:end, :), 'Color', cmatrix(1:1000:end, :));
                % pcshow(pcdepthtransformed, 'MarkerSize', 50);
                % hold on;
                pointclouds_depth{ind} = depthdata;
                % pcwrite(pcdepthtransformed, [folder, 'long_throw_depth/' , parsed{1}, '_transf.ply']);
                ind = ind + 1;
            end
        end
    end
end

% calculate mean rotation between HoloLens and Vicon coordinate systems
sumqs = sum(allqs, 2);
Qstar = sumqs / norm(sumqs);
Rstar = q2r(Qstar);

%% ---------- Plot results of depth pointcloud transforming  ---------- 

holdata = cell2mat(pointclouds_depth');

% crop and subsample pointclouds
minx = min(holdata(:,1));
miny = min(holdata(:,2));
minz = min(holdata(:,3));

maxx = max(holdata(:,1));
maxy = max(holdata(:,2));
maxz = max(holdata(:,3));

ptsH = [minx minx minx minx maxx maxx maxx maxx; 
       miny miny maxy maxy miny miny maxy maxy;
       minz maxz minz maxz minz maxz minz maxz]';

m = pcMatterport.Location;
indexes = logical(m(:, 1) >= -4) & logical(m(:, 1) <= 6);
ptsM = m(indexes,:);

pcptsM = pointCloud(ptsM(1:20:end, :));

indH =  logical(holdata(:, 1) >= -4) & logical(holdata(:, 1) <= 6);
ptsH = holdata(indH,:);
pcptsH = pointCloud(ptsH(1:20:end, :));


% plot pointclouds
figure();
pcshowpair(pcptsH, pcptsM, 'MarkerSize', 50);


% plot pointclouds with NN 
figure();
pcshowpair(pcptsH, pcptsM, 'MarkerSize', 50);
hold on;


% search NN of transformed HoloLens depth in Matterport
Mdl = KDTreeSearcher(pcptsM.Location);
[Idx, D] = knnsearch(Mdl, pcptsH.Location);

plot3(pcptsM.Location(Idx, 1), pcptsM.Location(Idx, 2),  pcptsM.Location(Idx, 3), 'or');
title('HoloLens and Matterport alignemnt depth PC  after tuning');

% plot errors 
figure();
plot(sort(D));
xlabel('Err id');
ylabel('Error [m]');
title('HoloLens and Matterport depth PC alignemnt errors after tuning');
%% ---------- Try depth pointcloud transforming using Vicon poses  ---------- 

count = 1;

for i = 1:length(files)
    % selecting correct files
    if ~startsWith(files(i).name, 'world_') && endsWith(files(i).name, '.ply') && ~endsWith(files(i).name,'_transf.ply')
        % parsing name and row number
        name = files(i).name;
        parsed = split(name, '.');
        poseID = ['long_throw_depth\' , parsed{1}, '.pgm'];
        row = depthposes(strcmp(depthposes.ImageFileName, poseID), :);
        % if file is in depth csv file
        if ~isempty(row)
            [~, rowindex] = ismember(poseID, depthposes.ImageFileName);
            
            % load HoloLesns depth pointcloud
            pcdepth = pcread([folder, 'long_throw_depth/', name]);
            depthdata = pcdepth.Location / -1000;
            pcdepth = pointCloud(depthdata);

            % get HoloLens pose matrices
            FrameToOrigin = [row.FrameToOrigin_m11 row.FrameToOrigin_m12 row.FrameToOrigin_m13 row.FrameToOrigin_m14;
                            row.FrameToOrigin_m21 row.FrameToOrigin_m22 row.FrameToOrigin_m23 row.FrameToOrigin_m24;
                            row.FrameToOrigin_m31 row.FrameToOrigin_m32 row.FrameToOrigin_m33 row.FrameToOrigin_m34;
                            row.FrameToOrigin_m41 row.FrameToOrigin_m42 row.FrameToOrigin_m43 row.FrameToOrigin_m44];

            CameraViewTransform = [row.CameraViewTransform_m11 row.CameraViewTransform_m12 row.CameraViewTransform_m13 row.CameraViewTransform_m14;
                                row.CameraViewTransform_m21 row.CameraViewTransform_m22 row.CameraViewTransform_m23 row.CameraViewTransform_m24;
                                row.CameraViewTransform_m31 row.CameraViewTransform_m32 row.CameraViewTransform_m33 row.CameraViewTransform_m34;
                                row.CameraViewTransform_m41 row.CameraViewTransform_m42 row.CameraViewTransform_m43 row.CameraViewTransform_m44];

            C2D = inv(CameraViewTransform)';
            D2C = CameraViewTransform';
            D2O = FrameToOrigin'; 
            O2D = inv(D2O);
            Rhi = D2O * C2D;
            oldRhi = Rhi(1:3, 1:3);
            mChw = oldRhi' * Rhi(1:3,4);
            Rv = euler2mat(Rvic(rowindex,:));
            Rhv_star = Rstar * Rv;
            corrections = hol2viccorr{6};
            
            
            %transform HoloLens pointcloud to Vicon using Vicon rotation
            scale = 1/optimizedPars.rho;
            rotation = optimizedPars.St * Rii';
            translation = T + Rv * optimizedPars.t(:,6) + corrections(rowindex,:)';
            depthdata = transformDepthDataToVicon(depthdata, Rhi, scale, rotation, translation);

            % save pointcloud
            cmatrix = ones(pcdepth.Count, 1) * cols{i};   
            pcDepthTransformed = pointCloud(depthdata(1:1000:end, :), 'Color', cmatrix(1:1000:end, :));
            pointclouds_depth_vicon{count} = depthdata;
            count = count + 1;

        end
    end
end

holdatavic = cell2mat(pointclouds_depth_vicon');

% crop and subsample HoloLens pointcloud
ptsH = holdatavic(indH,:);
pcptsH = pointCloud(ptsH(1:20:end, :));

% plot HoloLens and Matterport
figure();
pcshowpair(pcptsH, pcptsM, 'MarkerSize', 50);
hold on;

% find NN of HoloLens in Matterport
[Idx, Dvicon] = knnsearch(Mdl, pcptsH.Location);

plot3(pcptsM.Location(Idx, 1), pcptsM.Location(Idx, 2),  pcptsM.Location(Idx, 3), 'or');
title('HoloLens and Matterport alignemnt using Vicon');

% plot errors (NN distances) between HoloLens and Matterport pointclouds
figure();
plot(sort(Dvicon));
xlabel('Err id');
ylabel('Error [m]');
title('HoloLens and Matterport depth PC alignemnt errors after using Vicon');

%% -------------------------------------------------------------------------------------------

function depthdata = transformDepthDataToVicon(pointcloud, pointcloudTransform, scale, rotation, translation)
                tmp = pointcloudTransform * [pointcloud ones(size(pointcloud,1),1)]';
                tmp = (scale * rotation * tmp(1:3,:)) + translation;
                depthdata = tmp';
end

function drawCamera(C, R)
    cameraPointsInWorld = [[diag([1 1 -1]) * 0.05; [0 0 0]] ones(4,1)]; 
    p = zeros(4,3);
    for i = 1:4
        p(i,:) = (R * cameraPointsInWorld(i,1:3)' + C)';
    end
    hold on;
    plot3([p(4, 1), p(1, 1)], [p(4, 2), p(1, 2)], [p(4, 3), p(1, 3)], 'r', 'LineWidth', 1);
    plot3([p(4, 1), p(2, 1)], [p(4, 2), p(2, 2)], [p(4, 3), p(2, 3)], 'g', 'LineWidth', 1);
    plot3([p(4, 1), p(3, 1)], [p(4, 2), p(3, 2)], [p(4, 3), p(3, 3)], 'b', 'LineWidth', 1);
end


