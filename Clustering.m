function detectionClusters = clusterDetections(detections, vehicleSize)
% loop through every single detection point if list is from same sensor
N= numel(detections);
distances = zeros(N);
% and measure eucledian distance between all of them
for i = 1:N
    for j = i+1:N
        if detections{i}.SensorIndex ==  detections{j}.SensorIndex
            distances(i,j) = norm(detections{i}.Measurement(1:2) - detections{j}.Measurement(1:2));
        else 
            distances(i,j) = inf;
        end
    end
end

leftToCheck = 1:N;
i =0;
detectionClusters= cell(N:1);
while ~isempty(leftToCheck)
    % remove the detections that are in the same cluster as the one under
    % consideration
    underConsideration = leftToCheck(1);
    clusterInds = (distances(underConsideration, leftToCheck) < vehicleSize);
    detInds = leftToCheck(clusterInds);
    clusterDets = [detections{detInds}];
    clusterMeas = [clusterDets.Measurement];
    meas = mean(clusterMeas , 2);
    meas2D = [mean(1:2);meas(4:5)];
    i = i + 1;
    detectionClusters{i} = detections{detInds(1)};
    detectionClusters{i}.Measurement = meas2D;
    leftToCheck(clusterInds) = []; % removing cluster indexes that have already been detection or gone through
end
detectionClusters(i+1:end) = [];
%The while loop above consists of the following steps:

%Pick the first detection in the checklist and check for its clustering neighbors.
%If the distance between the first pick and remaining detections is less than the vehicle size, then group those detections and their respective radar sensor measurements, including range and velocity.
%For the group, take the mean of the range and velocity measurements.
%Note: the radar measurement vector has 6 values - where range and velocity for x and y coordinates reside at indices 1,2, 4, and 5: [x, y, -, Vx, Vy, -]

%Create a new Cluster ID. Then, assign all the group detections to the same ID.
%Further, assign cluster, the mean range, and velocity.
%In the end, delete from the list the detections which have already been assigned to a cluster.
%Keep repeating the process until the detection list is empty.