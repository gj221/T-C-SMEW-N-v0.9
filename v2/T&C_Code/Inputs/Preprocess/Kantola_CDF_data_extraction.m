% Create the main structure
data = struct();
data.version = [4, 2];

% Create the axes collection
data.axesColl = struct('name', 'XY', ...
                       'type', 'XYAxes', ...
                       'isLogX', false, ...
                       'isLogY', false, ...
                       'noRotation', false, ...
                       'calibrationPoints', []);

% Create and populate the Blueridge dataset
blueridge = struct('name', 'Blueridge', ...
                   'axesName', 'XY', ...
                   'colorRGB', [0, 0, 255, 255], ...
                   'metadataKeys', [], ...
                   'data', []);

blueridge.data = [
    struct('x', 85.6, 'y', 343.2, 'value', [0, 0]),
    struct('x', 85.6, 'y', 240.24, 'value', [250, 30]),
    struct('x', 171.2, 'y', 154.44, 'value', [500, 55]),
    struct('x', 342.4, 'y', 75.504, 'value', [1000, 78]),
    struct('x', 684.8, 'y', 27.456, 'value', [2000, 92]),
    struct('x', 1027.2, 'y', 6.864, 'value', [3000, 98]),
    struct('x', 1369.6, 'y', 0, 'value', [4000, 100])
];

data.datasetColl = {blueridge};

% Convert to JSON and save to file
jsonString = jsonencode(data, 'PrettyPrint', true);
fid = fopen('blueridge_CDF_data.json', 'w');
if fid == -1
    error('Cannot create JSON file');
end
fwrite(fid, jsonString, 'char');
fclose(fid);

disp('JSON file has been created: blueridge_data.json');