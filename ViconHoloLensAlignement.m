croppingStart = 0;
croppingEnd = 100;
Binit = 0;
BMax = 100;
%folder = './Vicon_session_2020_12_02/HoloLensRecording__2020_12_02__12_57_18/';
folder = './Vicon_session_2020_12_02/HoloLensRecording__2020_12_02__12_54_39/';
vicom = readtable('./Vicon_session_2020_12_02/hololens_seq03.txt');
pvhololens = readtable([folder, 'pv.csv']);
useAllCameras = true;

vlcll = readtable([folder, 'vlc_ll.csv']);
vlclf = readtable([folder, 'vlc_lf.csv']);
vlcrf = readtable([folder, 'vlc_rf.csv']);
vlcrr = readtable([folder, 'vlc_rr.csv']);
longthrowdepth = readtable([folder, 'long_throw_depth.csv']);

if useAllCameras
   allhololens = [pvhololens; vlcll; vlclf; vlcrf; vlcrr; longthrowdepth];
   allhololens = sortrows(pvhololens, 'Timestamp');
end


Xmin = 1.316;
Xmax = 1.445;
Ymin = -1.2105;
Ymax = -1.145;
Zmin = 1.212;
Zmax = 1.275;
%pri zmene sekvence se musi zmenit i folder dole!


% croppingStart = 100;
% croppingEnd = 0;
% Binit = 90;
% BMax = 150;
% vicom = readtable('./Vicom_2020_08_20/2020_08_20_Vicon_HoloLens_Session/seq1.txt');
% pvhololens = readtable('./Vicom_2020_08_20/HoloLensRecording__2020_08_20__08_36_27/pv.csv');
% folder = './Vicom_2020_08_20/HoloLensRecording__2020_08_20__08_36_27/';


indexes = vicom.Var4(:) ~= 1;

%% --------------- Odstraneni chyb
% for i = 1:size(vicom,1)
%     if vicom.Var5(i) > Xmin && vicom.Var5(i) < Xmax && vicom.Var6(i) > Ymin && vicom.Var6(i) < Ymax && vicom.Var7(i) > Zmin && vicom.Var7(i) < Zmax 
%         indexes(i) = 0;        
%     end
%     
%     
% end
%% --------------------------------
figure();

pcVicom = pointCloud([vicom.Var5(indexes), vicom.Var6(indexes), vicom.Var7(indexes)]);
ViconRot = [vicom.Var5(indexes), vicom.Var6(indexes), vicom.Var7(indexes)];
pcHoloLens = pointCloud([allhololens.Position_X(croppingStart+1:end-croppingEnd), allhololens.Position_Y(croppingStart+1:end-croppingEnd), allhololens.Position_Z(croppingStart+1:end-croppingEnd)]);
pcshowpair(pcHoloLens, pcVicom, 'MarkerSize', 50);

axis equal;
xlabel("x");
ylabel("y");
zlabel("z");
title('Hololens and Vicom tracking before registration');


%ax = 269;
besterr = realmax;
bestax = 0;
bestay = 0;
bestaz = 0;
bestrmse =0;
R = diag([1 1 1]);
res = cell(6,100);


MdlICP = KDTreeSearcher(pcVicom.Location);

for i = 1:100
    ay = interp1([0,1],[0,360],rand(1));
    az = interp1([0,1],[0,360],rand(1));
    ax = interp1([0,1],[0,360],rand(1));


    Rx = [1 0 0 0; 0 cos(ax) -sin(ax) 0; 0 sin(ax) cos(ax) 0; 0 0 0 1];
    Ry = [cos(ay) 0 sin(ay) 0; 0 1 0 0; -sin(ay) 0 cos(ay) 0; 0 0 0 1];
    Rz = [cos(az) -sin(az) 0 0; sin(az) cos(az) 0 0; 0 0 1 0; 0 0 0 1];

    tmpR = Rx * Ry * Rz;

    tform_rotate = affine3d(tmpR);
    ptrotHolol = pctransform(pcHoloLens,tform_rotate);

    [tform,tmphololensReg, rmse] = pcregrigid(ptrotHolol, pcVicom, 'MaxIterations', 100000, 'Tolerance', [0.001, 0.0009]);
    [Idx, D] = knnsearch(MdlICP, tmphololensReg.Location);
    %groups = kmeans(D, 2); %nemusi se vubec definovat, nebo mean(res) a mean(res) * 100
    %idxmin = groups(D == min(D));
    %err = sum(D(groups == idxmin));
    err = sum(D);
    
    res{1, i} = ax;
    res{2, i} = ay;
    res{3, i} = az;
    res{4, i} = err;
    res{5, i} = tform;
    res{6, i} = tmphololensReg;
    
    if err < besterr
        besterr = err;
        besttform = tform;
        hololensReg = tmphololensReg;
        bestax = ax;
        bestay = ay;
        bestaz = az;
        bestrmse = rmse;
        R = tmpR;
        bestD = D;
    end
% figure();
% pcshowpair(tmphololensReg, pcVicom, 'MarkerSize', 50);
% axis equal;
% xlabel("x");
% ylabel("y");
% zlabel("z");
% title(['Registered Hololens tracking and Vicom tracking, ax = ', num2str(ax), ', ay = ', num2str(ay), ', az = ', num2str(az), ' rmse = ', num2str(rmse)]);

end

% allerrs = [res{4,:}];
% allerrs = round(allerrs, 5);
% [GC, GR] = groupcounts(allerrs');
% err = GR(GC == max(GC));
% err_indexes = find(allerrs == err);
% 
% ax = res{1,err_indexes(1)};
% ay = res{2,err_indexes(1)};
% az = res{3,err_indexes(1)};
% tform = res{5,err_indexes(1)};
% hololensReg = res{6,err_indexes(1)};

% Rx = [1 0 0 0; 0 cos(ax) -sin(ax) 0; 0 sin(ax) cos(ax) 0; 0 0 0 1];
% Ry = [cos(ay) 0 sin(ay) 0; 0 1 0 0; -sin(ay) 0 cos(ay) 0; 0 0 0 1];
% Rz = [cos(az) -sin(az) 0 0; sin(az) cos(az) 0 0; 0 0 1 0; 0 0 0 1];

% R = Rx * Ry * Rz;

figure();
pcshowpair(hololensReg, pcVicom, 'MarkerSize', 50);
axis equal;
xlabel("x");
ylabel("y");
zlabel("z");
title(['Registered Hololens tracking and Vicom tracking, ax = ', num2str(ax), ', ay = ', num2str(ay), ', az = ', num2str(az)]);
tform_rotate = affine3d(R);
%% --------------------------------------------------------------

numHolol = size(pvhololens,1);
numVicom = size(vicom,1);

params = [];
vals = [];
initvals = [];

alldata = {pvhololens, vlcll, vlclf, vlcrf, vlcrr, longthrowdepth};

pointclouds = cell(6,1);
cs = cell(6,1);
hol = cell(6,1);
tform_rotate = affine3d(R);

for i = 1:6
    if i == 6
        pointclouds{i} = pointCloud([alldata{i}.Position_X(2:end), alldata{i}.Position_Y(2:end), alldata{i}.Position_Z(2:end)]);
        pointclouds{i} = pctransform(pctransform(pointclouds{i},tform_rotate), besttform);
        hol{i} = pointclouds{i}.Location;
        cs{i} = round((alldata{i}.Timestamp(2:end) - min(alldata{1}.Timestamp)) /10^5) + 1;
    else
        pointclouds{i} = pointCloud([alldata{i}.Position_X, alldata{i}.Position_Y, alldata{i}.Position_Z]);
        pointclouds{i} = pctransform(pctransform(pointclouds{i},tform_rotate), besttform);
        hol{i} = pointclouds{i}.Location;
        cs{i} = round((alldata{i}.Timestamp - min(alldata{1}.Timestamp)) /10^5) + 1;
    
    end
end
%hol = hololensReg.Location;
minB = 0;

% timediff = pvhololens.Timestamp(2 + croppingStart:end - croppingEnd) - pvhololens.Timestamp(1 + croppingStart:end - croppingEnd -1);
% timediff = timediff / 10^7;
% cs = cumsum([0; timediff]);
% cs  = round(cs * 100) + 1;

for B = Binit:BMax 
    vic = cell(6,1);
    Rvic = cell(6,1);
    vicsize = size(pcVicom.Location, 1);
    for k = 1:6
        i = cs{k} + B; 
        i(i > vicsize) = [];
        vic{k} = pcVicom.Location(i, :);
        Rvic{k} = ViconRot(i, :);
    end
    
    fun = @(par)closerTransformation([[par(1) par(2) par(3)]', [par(4) par(5) par(6)]', ...
    [par(7) par(8) par(9)]', [par(10) par(11) par(12)]', [par(13) par(14) par(15)]', ...
    [par(16) par(17) par(18)]'], par(19), euler2mat([par(20) par(21) par(22)]), [par(23) par(24) par(25)]', ...
    vic, Rvic, hol);
    initpar = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0];
    options = optimset( 'MaxFunEvals', 1000 * 25);
    %fun(initpar);
    [bestpar, bestval] = fminsearch(fun, initpar, options);
    
    %t, rho, S, d, C, R, D
    initvals = [initvals closerTransformation([[initpar(1) initpar(2) initpar(3)]', [initpar(4) initpar(5) initpar(6)]', ...
        [initpar(7) initpar(8) initpar(9)]', [initpar(10) initpar(11) initpar(12)]', [initpar(13) initpar(14) initpar(15)]', ...
        [initpar(16) initpar(17) initpar(18)]'], initpar(19), euler2mat([initpar(20) initpar(21) initpar(22)]), [initpar(23) initpar(24) initpar(25)]', vic, Rvic, hol)];
    
    params = [params bestpar];
    vals = [vals bestval];
    
    if B == Binit
      minpar = bestpar;
      minval = bestval;
      minB = B;
    end
    
    if bestval < minval
      minpar = bestpar;
      minval = bestval; 
      minB = B;
    end
    fprintf(['For B = ', num2str(B), ' the optimized value is ',  num2str(bestval),' parameters are:', '\n' ]);
    bestpar
    
    %t, rho, S, d, C, R, D
end

fprintf(['Best B is ', num2str(minB), ' the optimized value is ',  num2str(minval),' parameters are:', '\n' ]);
minpar

%% --------------------------------Commented because for all camera this doesn't make sense-----------------------------------------------------------
% %pcVicom = pointCloud(vic);
% holtrans = zeros(hololensReg.Count, 3);
% hol2vicon = zeros(hololensReg.Count, 3);
% j = uint64(cs + minB); 
% vic = pcVicom.Location(j, :);
% 
% Rvic = ViconRot(j, :);
% 
% for i = 1:hololensReg.Count
%     holtrans(i,:) = (euler2mat([minpar(5) minpar(6) minpar(7)])' * ((1/minpar(4)) * hol(i,:)') - [minpar(8) minpar(9) minpar(10)]')'; %+ Rvic * t cary jinou barvou;
%     hol2vicon(i, :) = (holtrans(i, :)' + euler2mat(Rvic(i, :)) * [minpar(1) minpar(2) minpar(3)]')';
% end
% 
% pcHol = pointCloud(holtrans);
% %% ---------------------------------------------------------------------------
% 
% 
% 
% 
% figure();
% %%
% 
% pcshowpair(pcHol, pcVicom, 'MarkerSize', 50);
% axis equal;
% xlabel("x");
% ylabel("y");
% zlabel("z");
% title('Registered Hololens tracking and Vicom tracking after tuning');
% 
% %%
% hold on;
% %kresleni car
% 
% cmatrix = ones(size(hol2vicon)).*[0 1 1];
% pcHol2Vicon = pointCloud(hol2vicon, 'Color', cmatrix);
% pcshow(pcHol2Vicon, 'MarkerSize', 50);
% hold on;
% 
% for i = 1:hololensReg.Count
%     plot3([vic(i, 1), holtrans(i,1)], [vic(i, 2), holtrans(i,2)], [vic(i, 3), holtrans(i,3)], 'w');
%     plot3([vic(i, 1), hol2vicon(i,1)], [vic(i, 2), hol2vicon(i,2)], [vic(i, 3), hol2vicon(i,3)], 'r');
% end
% 
% 
% legend('\color{white} Transformed HoloLens (without Rvic * t)','\color{white} Vicon', '\color{white} Transformed Hololens (with Rvic * t)', '\color{white} Error between HoloLens without Rvic * t and Vicon', '\color{white} Error between HoloLens with Rvic * t and Vicon');
% 
% 
% %%
% figure();
% plot(Binit:BMax , vals, 'b', Binit:BMax , initvals, 'r');
% 
% holCheck = pcHoloLens.Location;
% 
% rho = minpar(4);
% St = euler2mat([minpar(5) minpar(6) minpar(7)])';
% Rii = R(1:3, 1:3) * tform.Rotation;
% T = St * tform.Translation' - [minpar(8) minpar(9) minpar(10)]';
% 
% for i = 1:hololensReg.Count
%     holCheck(i,:) = 1/rho * St *  (holCheck(i,:) * Rii)' + T;
% end
% 
% pcHolCheck = pointCloud(holCheck);
% 
% figure();
% pcshowpair(pcHolCheck, pcVicom, 'MarkerSize', 50);
% axis equal;
% xlabel("x");
% ylabel("y");
% zlabel("z");
% title('Registered Hololens tracking and Vicom tracking after tuning');
% 
% 
% %%
% hold on;
% 
% plot3(vic(:, 1), vic(:, 2), vic(:, 3), 'ro');
% 
% %% ------------------- IN MATTERPORT
% pcMatterport = pcread('D:/Documents/Work/B-315/matterport2vicon.ply');
% figure();
% 
% pcshowpair(pcHolCheck, pcVicom, 'MarkerSize', 50);
% hold on;
% pcshow(pcMatterport, 'MarkerSize', 50);
% axis equal;
% xlabel("x");
% ylabel("y");
% zlabel("z");
% title('Transformed HoloLens and Vicon data in transformed Matterport data');
% legend('\color{white} Transformed HoloLens','\color{white} Vicon', '\color{white} Matterport');
% 
% %% ------------------------ graphs/histograms
% 
% errorHolTrans = sum((holtrans - vic).^2, 2);
% errorHol2Vicon = sum((hol2vicon - vic).^2, 2);
% 
% figure()
% hold on;
% grid on;
% xlabel("pose id");
% ylabel("error (m)");
% 
% plot(errorHolTrans);
% plot(errorHol2Vicon);
% title('Graph of errors');
% legend('Transformed HoloLens (without Rvic * t)','Transformed HoloLens (with Rvic * t)');
% 
% 
% figure()
% hold on;
% grid on;
% xlabel("error (m)");
% ylabel("number of occurence");
% histogram(errorHolTrans);
% histogram(errorHol2Vicon);
% 
% title('Histogram of errors');
% legend('Transformed HoloLens (without Rvic * t)','Transformed HoloLens (with Rvic * t)');
%% ------- Grafy
figure();
plot(Binit:BMax , vals, 'b', Binit:BMax , initvals, 'r');
legend('Optimized function value', 'Initial function value');

St = euler2mat([minpar(20) minpar(21) minpar(22)])';
rho = minpar(19);
d = [minpar(23) minpar(24) minpar(25)]';
t = [[minpar(1) minpar(2) minpar(3)]' ...
[minpar(4) minpar(5) minpar(6)]' ...
[minpar(7) minpar(8) minpar(9)]' ...
[minpar(10) minpar(11) minpar(12)]' ...
[minpar(13) minpar(14) minpar(15)]' ...
[minpar(16) minpar(17) minpar(18)]'];

for m = 1:6

hol2vicon = zeros(size(hol{m},1), 3);

j = cs{m} + minB; 
j(j > pcVicom.Count) = [];
vic = pcVicom.Location(j, :);
victrans = zeros(size(hol{m},1), 3);
 
Rvic = ViconRot(j, :);

for i = 1:size(hol{m},1)
    holtrans = St * ((1/rho) * hol{m}(i,:)') + d; 
    hol2vicon(i, :) = (holtrans)';
    victrans(i, :) = (vic(i, :)' + euler2mat(Rvic(i, :)) * t(:,m))';
end

pcHol = pointCloud(hol2vicon);
 
figure();

pcshowpair(pcHol, pcVicom, 'MarkerSize', 50);
axis equal;
xlabel("x");
ylabel("y");
zlabel("z");
title(['Registered Hololens tracking and Vicom tracking after tuning for ', num2str(m), 'th camera']);



figure();
axis equal;
xlabel("x");
ylabel("y");
zlabel("z");
title(['Registered Hololens tracking and Vicom tracking after tuning for ', num2str(m), 'th camera']);
pcshowpair(pcHol, pcVicom, 'MarkerSize', 50);
hold on;

for i = 1:pcHol.Count
    plot3([victrans(i, 1), hol2vicon(i,1)], [victrans(i, 2), hol2vicon(i,2)], [victrans(i, 3), hol2vicon(i,3)], 'r');
end
hold on;
pcshow(pointCloud(hol{m}), 'MarkerSize', 50);
pcshow(pointCloud(victrans, 'Color', repmat([0 1 1], length(victrans), 1)),'MarkerSize', 50);
legend('\color{white} Transformed HoloLens','\color{white} Vicon', '\color{white} Error between HoloLens' , '\color{white} Hololens before tuning', '\color{white} Shifted Vicon');

dist = vecnorm((victrans - hol2vicon)');
distgroups = kmeans(dist', 2);
idxmin = distgroups(dist == min(dist));
meanerr = mean(dist(distgroups == idxmin));
fprintf(['For ', num2str(m), 'th camera the mean distance is: ',  num2str(meanerr), '\n' ]);
hol2viccorr{m} = victrans - hol2vicon;

end



%% ------------               -------- Depth data

files=dir([folder, 'long_throw_depth/']);
depthposes = readtable([folder, 'long_throw_depth.csv']);
depthposes(1,:) = [];
% figure();
% hold on;
% grid on;
% xlabel("x");
% ylabel("y");
% zlabel("z");
% title('HoloLens Depth Data pointclouds');
pcMatterport = pcread('D:/Documents/Work/B-315/matterport2vicon.ply');
% pcshow(pcMatterport, 'MarkerSize', 50);
pointclouds_depth = {};

ind = 1;
cols = getColors(length(files));

j = cs{6} + minB; 
Rvic = ViconRot(j, :);
vic = pcVicom.Location(j, :);

Rii = R(1:3, 1:3) * besttform.Rotation;
T = St * besttform.Translation' - d; % mozna -d mozna +d
qi = [0 0 0 0];


%     holtrans(i,:) = (euler2mat([minpar(5) minpar(6) minpar(7)])' * ((1/minpar(4)) * hol(i,:)') - [minpar(8) minpar(9) minpar(10)]')'; %+ Rvic * t cary jinou barvou;
%     hol2vicon(i, :) = (holtrans(i, :)' + euler2mat(Rvic(i, :)) * [minpar(1) minpar(2) minpar(3)]')';

for i = 1:length(files)
   if ~startsWith(files(i).name, 'world_') && endsWith(files(i).name, '.ply') && ~endsWith(files(i).name,'_transf.ply')
%        if ind == 10
       name = files(i).name;
       parsed = split(name, '.');
       
       
       poseID = ['long_throw_depth\' , parsed{1}, '.pgm'];
       row = depthposes(strcmp(depthposes.ImageFileName, poseID), :);
       if ~isempty(row)
       [~, rowindex] = ismember(poseID, depthposes.ImageFileName);
       cameraT = [row.Position_X row.Position_Y row.Position_Z];
       cameraQ = [row.Orientation_W row.Orientation_X row.Orientation_Y row.Orientation_Z];
       FrameToOrigin = [row.FrameToOrigin_m11 row.FrameToOrigin_m12 row.FrameToOrigin_m13 row.FrameToOrigin_m14;
                        row.FrameToOrigin_m21 row.FrameToOrigin_m22 row.FrameToOrigin_m23 row.FrameToOrigin_m24;
                        row.FrameToOrigin_m31 row.FrameToOrigin_m32 row.FrameToOrigin_m33 row.FrameToOrigin_m34;
                        row.FrameToOrigin_m41 row.FrameToOrigin_m42 row.FrameToOrigin_m43 row.FrameToOrigin_m44];
       
       CameraViewTransform = [row.CameraViewTransform_m11 row.CameraViewTransform_m12 row.CameraViewTransform_m13 row.CameraViewTransform_m14;
                              row.CameraViewTransform_m21 row.CameraViewTransform_m22 row.CameraViewTransform_m23 row.CameraViewTransform_m24;
                              row.CameraViewTransform_m31 row.CameraViewTransform_m32 row.CameraViewTransform_m33 row.CameraViewTransform_m34;
                              row.CameraViewTransform_m41 row.CameraViewTransform_m42 row.CameraViewTransform_m43 row.CameraViewTransform_m44];
      
%        cameraPoints = [[diag([1 1 1]) * 0.3; [0 0 0]] ones(4,1)];              
%        transformedCamPoints = zeros(5,3);
%        Rit = [[Rii'; 0 0 0] [0; 0; 0; 1]];
%        CamTran = [[inv(CameraViewTransform(1:3, 1:3));0 0 0] [0; 0; 0; 1]];
%        F2O = [[FrameToOrigin(1:3, 1:3);0 0 0] [0; 0; 0; 1]];
%        
%        D2C = inv(CameraViewTransform)';
%        C2D = inv(D2C);
%        translation = [FrameToOrigin(4,1) -FrameToOrigin(4,2) -FrameToOrigin(4,3) FrameToOrigin(4,4)]; %+ CameraViewTransform(4,:);

        C2D = inv(CameraViewTransform)';
        D2C = inv(C2D);
        D2O = FrameToOrigin'; 
        O2D = inv(D2O);
        Rhi = D2O(1:3, 1:3) * C2D(1:3, 1:3);
        
%        for m = 1:4
%           cameraPoints(m, :) =  (D2C * cameraPoints(m, :)')' + translation;
%        end
%        drawCamera(transformedCamPoints);
  
       %text(transformedCamPoints(1,1), transformedCamPoints(1,2), transformedCamPoints(1,3), parsed{1});
       pcdepth = pcread([folder, 'long_throw_depth/', name]);
       depthdata = pcdepth.Location / -1000;
       pcdepth = pointCloud(depthdata);
       corrections = hol2viccorr{6};
       %pcdepth = pctransform(pctransform(pcdepth,tform_rotate), besttform); % timhle se zbavime Rii
       %cmatrix = ones(size(depthdata)) .* cols{i};
       %cmatrix = ones(size(depthdata)) .* [0 1 0];
       %Rx = [1 0 0 0; 0 cos(pi) -sin(pi) 0; 0 sin(pi) cos(pi) 0; 0 0 0 1];
       %Rx3 = [1 0 0; 0 cos(pi) -sin(pi); 0 sin(pi) cos(pi)];
       if(norm(corrections(rowindex,:)) < 0.2)
           Rhv = St * Rii' * Rhi;
           Rv = euler2mat(Rvic(rowindex,:));
           qdelta = r2q(Rhv * Rv'); 
           qi = qi + qdelta;

           Rdelta = q2r(qdelta);
           % chceme, aby platilo Rhv = Rdelta * Rv;
           
           for k = 1:pcdepth.Count
               % depthdata to hololens world coordinates
               tmp = ((D2O) * (C2D) * [depthdata(k,:) 1]')'; 

               % hololens world -> 
               tmp = (1/rho * St * Rii' * tmp(1:3)') + T + euler2mat(Rv) * t(:,6) + corrections(rowindex,:)';
               %tmp = 1/rho * St * ((Rhi * depthdata(k,:)')' *  Rii)' + victrans(rowindex,:)';
               %tmp = (1/rho * St * Rii' * tmp(1:3)') + T + euler2mat(Rvic(rowindex,:)) * t(:,6);


               %tmp = euler2mat(Rvic(ind,:)) * depthdata(k,:)' + euler2mat(Rvic(ind,:)) * t(:,6) + vic(ind,:);
               depthdata(k, :) = [tmp(1) tmp(2) tmp(3)];
           end
        %        pcdepth1 = pctransform(pointCloud(depthdata), tform_rotate);
        %        pcdepth2 = pctransform(pcdepth1,tform);
        %        xyz = pcdepth2.Location;
        %        for k = 1:pcdepth2.Count
        %            xyz(k, :) = 
        %        end
               cmatrix = ones(pcdepth.Count, 1) * cols{i};   
               pcdepthtransformed = pointCloud(depthdata(1:1000:end, :), 'Color', cmatrix(1:1000:end, :));
        %        pcshow(pcdepthtransformed, 'MarkerSize', 50);
               hold on;
               pointclouds_depth{ind} = depthdata;
               %pcwrite(pcdepthtransformed, [folder, 'long_throw_depth/' , parsed{1}, '_transf.ply']);
        %        break;

        %        end
               ind = ind + 1;
           end
       end
   end
end

%%

%Qstar = qi / (ind - 1);
Qstar = qi / norm(qi);
Rstar = q2r(Qstar);




holdata = cell2mat(pointclouds_depth');
figure();
pcshow(pointCloud(holdata));
title('HoloLens transformed to Vicon')

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

pcshowpair(pcptsH, pcptsM, 'MarkerSize', 50);

Mdl = KDTreeSearcher(pcptsM.Location)
% [tform,hololFixed] = pcregrigid(pcptsH, pcptsM, 'MaxIterations', 10000, 'Verbose', true, 'Tolerance', [0.0001, 0.00009]);
% figure();
% pcshowpair(hololFixed, pcptsM, 'MarkerSize', 50);
% hold on;
% 
% ;
% [IdxICP, DICP] = knnsearch(Mdl, hololFixed.Location);
% 
% plot3(pcptsM.Location(IdxICP, 1), pcptsM.Location(IdxICP, 2),  pcptsM.Location(IdxICP, 3), 'or');
% title('HoloLens and Matterport alignemnt after ICP');
% 
% figure();
% plot(D);
% hold on;
% plot(DICP);
% title('Distances tuned vs tuned and ICP');
% legend('D', 'D ICP');
% 
% figure();
% histogram(D,100);
% hold on;
% histogram(DICP,100);
% title('Distances tuned vs tuned and ICP');
% legend('D', 'D ICP');
% 
% 
figure();
pcshowpair(pcptsH, pcptsM, 'MarkerSize', 50);
hold on;

[Idx, D] = knnsearch(Mdl, pcptsH.Location);

plot3(pcptsM.Location(Idx, 1), pcptsM.Location(Idx, 2),  pcptsM.Location(Idx, 3), 'or');
title('HoloLens and Matterport alignemnt without ICP');

%%

% %% -------------------------------- Porovnani timestampu -------------------------------------
% pvTimestamps = pvhololens.Timestamp(croppingStart+1:end-croppingEnd);
% viconTimestamps = vicom.Var1(indexes);
% viconMatchedTimestamps = viconTimestamps(j);
% depthTimestamps = depthposes.Timestamp;
% 
% figure();
% plot(pvTimestamps, '.');
% grid on;
% title("PV timestamps");
% 
% 
% figure();
% plot(viconTimestamps, '.');
% grid on;
% title("Original Vicon timestamps");
% 
% figure();
% plot(viconMatchedTimestamps, '.');
% grid on;
% title("Vicon timestamps matched to PV");
% 
% figure();
% plot(depthTimestamps, '.');
% grid on;
% title("Long throw depth timestamps");
% 
% 
% figure();
% plot(pvTimestamps, viconMatchedTimestamps);
% grid on;
% title("Mapping of Vicon timestamp to PV timestamps");
% xlabel("PV timestamps");
% ylabel("Mapped Vicon Timestamps");


count = 1;

for i = 1:length(files)
   if ~startsWith(files(i).name, 'world_') && endsWith(files(i).name, '.ply') && ~endsWith(files(i).name,'_transf.ply')
%        if ind == 10
       name = files(i).name;
       parsed = split(name, '.');
       
       
       poseID = ['long_throw_depth\' , parsed{1}, '.pgm'];
       row = depthposes(strcmp(depthposes.ImageFileName, poseID), :);
       if ~isempty(row)
       [~, rowindex] = ismember(poseID, depthposes.ImageFileName);

       pcdepth = pcread([folder, 'long_throw_depth/', name]);
       depthdata = pcdepth.Location / -1000;
       pcdepth = pointCloud(depthdata);

       FrameToOrigin = [row.FrameToOrigin_m11 row.FrameToOrigin_m12 row.FrameToOrigin_m13 row.FrameToOrigin_m14;
                row.FrameToOrigin_m21 row.FrameToOrigin_m22 row.FrameToOrigin_m23 row.FrameToOrigin_m24;
                row.FrameToOrigin_m31 row.FrameToOrigin_m32 row.FrameToOrigin_m33 row.FrameToOrigin_m34;
                row.FrameToOrigin_m41 row.FrameToOrigin_m42 row.FrameToOrigin_m43 row.FrameToOrigin_m44];
       
       CameraViewTransform = [row.CameraViewTransform_m11 row.CameraViewTransform_m12 row.CameraViewTransform_m13 row.CameraViewTransform_m14;
                              row.CameraViewTransform_m21 row.CameraViewTransform_m22 row.CameraViewTransform_m23 row.CameraViewTransform_m24;
                              row.CameraViewTransform_m31 row.CameraViewTransform_m32 row.CameraViewTransform_m33 row.CameraViewTransform_m34;
                              row.CameraViewTransform_m41 row.CameraViewTransform_m42 row.CameraViewTransform_m43 row.CameraViewTransform_m44];
      
        C2D = inv(CameraViewTransform)';
        D2C = inv(C2D);
        D2O = FrameToOrigin'; 
        O2D = inv(D2O);
        Rhi = D2O * C2D;
        oldRhi = Rhi(1:3, 1:3);
        mChw = oldRhi' * Rhi(1:3,4);
        Rhi = [Rstar * euler2mat(Rvic(rowindex,:)) Rstar * euler2mat(Rvic(rowindex,:)) * mChw; 0 0 0 1];
        
        corrections = hol2viccorr{6};
        
       for k = 1:pcdepth.Count
                          % depthdata to hololens world coordinates
               tmp = (Rhi * [depthdata(k,:) 1]')'; 
               %tmp = Rii * St' * tmp(1:3)';

               % hololens world -> 
               %tmp = (1/rho * St *  (tmp(1:3) * Rii)') + T + euler2mat(Rvic(rowindex,:)) * t(:,6) + corrections(rowindex,:)';

           tmp =  1/rho * tmp(1:3) + T + euler2mat(Rvic(rowindex,:)) * t(:,6) + corrections(rowindex,:)';
           depthdata(k, :) = [tmp(1) tmp(2) tmp(3)];
       end
%        pcdepth1 = pctransform(pointCloud(depthdata), tform_rotate);
%        pcdepth2 = pctransform(pcdepth1,tform);
%        xyz = pcdepth2.Location;
%        for k = 1:pcdepth2.Count
%            xyz(k, :) = 
%        end
       cmatrix = ones(pcdepth.Count, 1) * cols{i};   
       pcdepthtransformed = pointCloud(depthdata(1:1000:end, :), 'Color', cmatrix(1:1000:end, :));
%        pcshow(pcdepthtransformed, 'MarkerSize', 50);
       hold on;
       pointclouds_depth_vicon{count} = depthdata;
       count = count + 1;

       end
   end
end

holdatavic = cell2mat(pointclouds_depth_vicon');
figure();
pcshow(pointCloud(holdatavic));
title('HoloLens transformed to Vicon')


ptsH = holdatavic(indH,:);
pcptsH = pointCloud(ptsH(1:20:end, :));

figure();
pcshowpair(pcptsH, pcptsM, 'MarkerSize', 50);
hold on;

[Idx, D] = knnsearch(Mdl, pcptsH.Location);

plot3(pcptsM.Location(Idx, 1), pcptsM.Location(Idx, 2),  pcptsM.Location(Idx, 3), 'or');
title('HoloLens and Matterport alignemnt using Vicon');
%% -------------------------------------------------------------------------------------------
function res = closerTransformation(t, rho, S, d, C, Rv, D)
%     C-known
%     R(e)-known
%     t-unknown
%     rho-unknown
%     S(e) - unknown
%     D - known
%     d - unknown

    res = cell(6,1);

     for k = 1:6
%     for k = 6
        
        num = size(Rv{k},1);
        res{k} = zeros(num, 1);

        for i = 1:num
            R = euler2mat(Rv{k}(i,:));
            f = C{k}(i,:)' + R * t(:,k) - (1/rho) * S' * D{k}(i,:)' - d;
            res{k}(i) = norm(f)^2;
            %res = res + log(norm(f) + 1);
        end
    end
    res = cell2mat(res);
    
    idx = kmeans(res, 2, 'Start', [min(res); max(res)]); %nemusi se vubec definovat, nebo mean(res) a mean(res) * 100
    idxmin = idx(res == min(res));
%    idxmax = max(idx);
%     
%     figure();
%     plot(res(idx == idxmin), 'r.');
%     hold on;
%     plot(res(idx == idxmax), 'b.');
%     
    res = sum(res(idx == idxmin(1))) + 10 * (norm(t) + norm(d)); %pak fixnout t(:,k) na t
    

    
    %res = sum(res(idx == idxmin))  + 10 * log(norm(t) + 1) + 10 * log(norm(d) + 1);
    %res = sum(norm(C + R*t - 1/ro * S' * D - d));
end

function R = euler2mat(e)
    x = e(1);
    y = e(2);
    z = e(3);
    R = [(cos(y) * cos(z)) (-1 * cos(y) * sin(z)) sin(y);
         (cos(x) * sin(z) + sin(x) * sin(y) * cos(z)) (cos(x) * cos(z) - sin(x) * sin(y) * sin(z)) (-1 * sin(x) * cos(y));
         (sin(x) * sin(z) - cos(x) * sin(y) * cos(z)) (sin(x) * cos(z) + cos(x) * sin(y) * sin(z)) (cos(x) * cos(y))];

end

function drawCamera(p)
hold on;
plot3([p(4, 1), p(1, 1)], [p(4, 2), p(1, 2)], [p(4, 3), p(1, 3)], 'r', 'LineWidth', 1);
plot3([p(4, 1), p(2, 1)], [p(4, 2), p(2, 2)], [p(4, 3), p(2, 3)], 'g', 'LineWidth', 1);
plot3([p(4, 1), p(3, 1)], [p(4, 2), p(3, 2)], [p(4, 3), p(3, 3)], 'b', 'LineWidth', 1);
end


