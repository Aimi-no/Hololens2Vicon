function allPointcloudData = createPointclouds(alldata, numOfCameras, tform_rotate, besttform)
% transforms all HoloLens data, creates pointclouds and calculates timestamp steps(cs) used to tuning alignment
% as a reference timestamp ('the first HoloLens timestamp') we are using
% the first timestamp of the PV camera

    allPointcloudData.pointclouds = cell(numOfCameras,1);
    allPointcloudData.cs = cell(numOfCameras,1);
    allPointcloudData.hol = cell(numOfCameras,1);

    for i = 1:numOfCameras
        if i == numOfCameras
            allPointcloudData.pointclouds{i} = pointCloud([alldata{i}.Position_X(2:end), alldata{i}.Position_Y(2:end), alldata{i}.Position_Z(2:end)]);
            allPointcloudData.pointclouds{i} = pctransform(pctransform(allPointcloudData.pointclouds{i},tform_rotate), besttform);
            allPointcloudData.hol{i} = allPointcloudData.pointclouds{i}.Location;
            allPointcloudData.cs{i} = round((alldata{i}.Timestamp(2:end) - min(alldata{1}.Timestamp)) /10^5) + 1;
        else
            allPointcloudData.pointclouds{i} = pointCloud([alldata{i}.Position_X, alldata{i}.Position_Y, alldata{i}.Position_Z]);
            allPointcloudData.pointclouds{i} = pctransform(pctransform(allPointcloudData.pointclouds{i},tform_rotate), besttform);
            allPointcloudData.hol{i} = allPointcloudData.pointclouds{i}.Location;
            allPointcloudData.cs{i} = round((alldata{i}.Timestamp - min(alldata{1}.Timestamp)) /10^5) + 1;

        end
    end


end