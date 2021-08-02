croppingStart = 0;
croppingEnd = 100;
Binit = 0;
BMax = 100;

[vicom, pvhololens, alldata, allhololens] = loadHololensData('./Vicon_session_2020_12_02/HoloLensRecording__2020_12_02__12_54_39/');
indexes = vicom.Var4(:) ~= 1;

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
legend('\color{white} Hololens', '\color{white} Vicon');


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
        besttform_rotate = tform_rotate;
        R = tmpR;
        bestD = D;
    end
end

figure();
pcshowpair(hololensReg, pcVicom, 'MarkerSize', 50);
axis equal;
xlabel("x");
ylabel("y");
zlabel("z");
title(['Registered Hololens tracking and Vicom tracking, ax = ', num2str(ax), ', ay = ', num2str(ay), ', az = ', num2str(az)]);
legend('\color{white} Hololens', '\color{white} Vicon');

%% --------------------------------------------------------------

numHolol = size(pvhololens,1);
numVicom = size(vicom,1);

params = [];
vals = [];
initvals = [];



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

minB = 0;

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
errs = [];

for m = 1:6

hol2vicon = zeros(size(hol{m},1), 3);

j = cs{m} + minB; 
j(j > pcVicom.Count) = [];
vic = pcVicom.Location(j, :);
victrans = zeros(size(hol{m},1), 3);
 
Rvic = ViconRot(j, :);

for i = 1:size(hol{m},1)
    holtrans = St * ((1/rho) * hol{m}(i,:)') - d; 
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

pcshow(pointCloud(hol{m}, 'Color', repmat([1 1 0], length(hol{m}),1)), 'MarkerSize', 50);
pcshow(pointCloud(victrans, 'Color', repmat([0 1 1], length(victrans), 1)),'MarkerSize', 50);
hold on;

for i = 1:pcHol.Count
    if(norm(victrans(i,:) - hol2vicon(i,:)) < 0.2)
        errs = [errs norm(victrans(i,:) - hol2vicon(i,:))];
        plot3([victrans(i, 1), hol2vicon(i,1)], [victrans(i, 2), hol2vicon(i,2)], [victrans(i, 3), hol2vicon(i,3)], 'r');
    end
end

legend('\color{white} Transformed HoloLens','\color{white} Vicon', '\color{white} Hololens before tuning', '\color{white} Shifted Vicon', '\color{white} Error between HoloLens');
title(['Registered Hololens tracking and Vicom tracking after tuning for ', num2str(m), 'th camera']);

dist = vecnorm((victrans - hol2vicon)');
distgroups = kmeans(dist', 2);
idxmin = distgroups(dist == min(dist));
meanerr = mean(dist(distgroups == idxmin));
fprintf(['For ', num2str(m), 'th camera the mean distance is: ',  num2str(meanerr), '\n' ]);
hol2viccorr{m} = victrans - hol2vicon;

end

figure();
plot(sort(errs));
xlabel("err id");
ylabel("error[m]");
title('All errors between aligned HoloLens and Vicon');

%% ------------               -------- Depth data

files=dir([folder, 'long_throw_depth/']);
depthposes = readtable([folder, 'long_throw_depth.csv']);
depthposes(1,:) = [];

pcMatterport = pcread('D:/Documents/Work/B-315/matterport2vicon.ply');
pointclouds_depth = {};

ind = 1;
cols = getColors(length(files));

j = cs{6} + minB; 
Rvic = ViconRot(j, :);
vic = pcVicom.Location(j, :);

Rii = R(1:3, 1:3) * besttform.Rotation;
T = St * besttform.Translation' - d; % mozna -d mozna +d
qi = [0 0 0 0];
allqs = [];
figure();

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
      

        C2D = inv(CameraViewTransform)';
        D2C = inv(C2D);
        D2O = FrameToOrigin'; 
        O2D = inv(D2O);
        Rhi = D2O * C2D;
        
        oldRhi = Rhi(1:3, 1:3);
        mChw = oldRhi' * Rhi(1:3,4);

       pcdepth = pcread([folder, 'long_throw_depth/', name]);
       depthdata = pcdepth.Location / -1000;
       pcdepth = pointCloud(depthdata);
       corrections = hol2viccorr{6};
       %pcdepth = pctransform(pctransform(pcdepth,tform_rotate), besttform); % timhle se zbavime Rii
       if(norm(corrections(rowindex,:)) < 0.2)
           Rv2Rvm = euler2mat(Rvic(rowindex,:));
           Rhv = St * Rii' * oldRhi;
           Rdelta = Rhv * Rv';
           qdelta = r2q(Rdelta); 
           allqs = [allqs qdelta];
           allRs{ind} =  Rdelta;
           qi = qi + qdelta;

           Rhv_star = Rdelta * Rv';
           
           Rstar_test = Rv' * St * Rii' * oldRhi';
           Rstarq{ind} = r2q(Rstar_test);
           % chceme, aby platilo Rhv = Rdelta * Rv';
%            drawCamera(-mChw, Rhv_star);
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
figure();
for i = 1:length(Rvic)
    drawCamera(vic(i,:)', euler2mat(Rvic(i,:)))
end
axis equal;
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

Mdl = KDTreeSearcher(pcptsM.Location);
% %
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
% % 
% % 
figure();
pcshowpair(pcptsH, pcptsM, 'MarkerSize', 50);
hold on;

[Idx, D] = knnsearch(Mdl, pcptsH.Location);

plot3(pcptsM.Location(Idx, 1), pcptsM.Location(Idx, 2),  pcptsM.Location(Idx, 3), 'or');
title('HoloLens and Matterport alignemnt depth PC  after tuning');

figure();
plot(sort(D));
xlabel('Err id');
ylabel('Error [m]');
title('HoloLens and Matterport depth PC alignemnt errors after tuning');
%%

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
        Rv = euler2mat(Rvic(rowindex,:));
        
        %Rhi = [Rstar * euler2mat(Rvic(rowindex,:)) Rstar * euler2mat(Rvic(rowindex,:)) * mChw; 0 0 0 1];
        Rhv_star = Rstar * Rv;
        corrections = hol2viccorr{6};
        
       for k = 1:pcdepth.Count
                          % depthdata to hololens world coordinates
               tmp = (Rhi * [depthdata(k,:) 1]')';

               % hololens world -> vicon
               tmp = (1/rho * St *  (tmp(1:3) * Rii)') + T + euler2mat(Rvic(rowindex,:)) * t(:,6) + corrections(rowindex,:)';
           %tmp = 1/rho * Rhv_star * (depthdata(k,:)' + mChw) + T + euler2mat(Rv) * t(:,6) + corrections(rowindex,:)'
           %tmp =  1/rho * tmp(1:3) + T + euler2mat(Rvic(rowindex,:)) * t(:,6) + corrections(rowindex,:)';
           depthdata(k, :) = [tmp(1) tmp(2) tmp(3)];
       end

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
function [vicom, pvhololens, alldata, allhololens] = loadHololensData(folder)
vicom = readtable('./Vicon_session_2020_12_02/hololens_seq03.txt');
pvhololens = readtable([folder, 'pv.csv']);
useAllCameras = true;

vlcll = readtable([folder, 'vlc_ll.csv']);
vlclf = readtable([folder, 'vlc_lf.csv']);
vlcrf = readtable([folder, 'vlc_rf.csv']);
vlcrr = readtable([folder, 'vlc_rr.csv']);
longthrowdepth = readtable([folder, 'long_throw_depth.csv']);
allhololens = []'

if useAllCameras
   allhololens = [pvhololens; vlcll; vlclf; vlcrf; vlcrr; longthrowdepth];
   allhololens = sortrows(pvhololens, 'Timestamp');
end

alldata = {pvhololens, vlcll, vlclf, vlcrf, vlcrr, longthrowdepth};

end


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


