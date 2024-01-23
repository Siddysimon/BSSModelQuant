%headNames = ["TIME041", "TIME045", "TIME052", "TIME057", "TIME072", "TIME081", "TIME082", "TIME085"];
% headNames = ["TIME052", "TIME072", "TIME081", "TIME082", "TIME085"];
headNames = ["TIME081"];


EAngles=[];
EMags = [];
ROIs = [];
Segs = [];
for i =1:length(headNames)
    %timeI = strrep(time20I, "TIME020", headNames(i));
    %timeF = strrep(time20F, "TIME020", headNames(i));
    configPath = strrep("/Volumes/spiral/BSSLab/MEQ/TIME020/setup/configs.mat", "TIME020", headNames(i));
    configs = load(configPath).configs;
    try
        electrodes = configs{1}(1,:);
    catch
        warning("config load didn't work");
    end 

    pre = strrep(preT20, "TIME020", headNames(i));
    preFit = strrep(preT20Fit, "TIME020", headNames(i));
    post = strrep(postT20, "TIME020", headNames(i));
    geo = strrep(geoT20Path, "TIME020", headNames(i));
    geoPre = strrep(geoT20PrePath, "TIME020", headNames(i));
   [EAngle, EMagnitude, ROI, Segmentation] = Analyze(electrodes, electrodes, geo, preFit, post, timeI, timeF, ROI_radii, headNames(i));
    EAngles = [EAngles; EAngle];
    EMags = [EMags; EMagnitude];
    ROIs = [ROIs; ROI];
    Segs = [Segs; Segmentation];
    EPost = load(post).E_vector;
    EPre = load(pre).E_vector;
    EPreFit = load(preFit).E_vector;
    GeometryPost = load(geo);
    GeometryPre = load(geoPre);

    brain_blocking(headNames(i), GeometryPost, GeometryPre, ROI, configs, EPost, EPreFit, EPre);
end
