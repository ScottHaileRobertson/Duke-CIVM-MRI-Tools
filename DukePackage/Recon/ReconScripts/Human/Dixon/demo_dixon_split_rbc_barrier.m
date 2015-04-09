%% Load data
disp('Load Ventilation volume');
vent_path = filepath('C:\Users\Scott\Desktop\subject002_065\P34304_gas_recon.mat')
load(vent_path);

disp('Load Dissolved volume');
dissolved_path = filepath('C:\Users\Scott\Desktop\subject002_065\P34304_dissolved_recon.mat')
load(dissolved_path);

disp('Load Proton volume');
proton_path = filepath('C:\Users\Scott\Desktop\subject002_065\P32256_bhute_recon.mat')
load(proton_path);

%% Find spectroscopic ratio
spec_ratio = 

%% Create segmentation mask
nClusters = 3;
dims = size(uteVol);
disp('Clustering');
[cluster_idx, cluster_center] = kmeans(abs(uteVol(:)),nClusters,'distance','sqEuclidean', ...
    'Replicates',3);
cluster_idx = reshape(cluster_idx,dims);

disp('Applying morphologic operations');


disp('Finding connected components');
[CC, NUM] = bwlabeln(cluster_idx==1, 18);
imslice(abs(uteVol),'UTE');
imslice(CC,'Cluster Indices');
CC_lung = input('Enter index of lung: ');

% Create mask
lung_mask = (CC == CC_lung);

%% Use the inverse cotangent to find where the ratio of real to imaginary 
% is the correct ratio
nPhases = 365;
phaseVec = linspace(0,2*pi,nPhases+1);
phaseVec = phaseVec(1:(end-1));
ratios = zeros(1,nPhases);
ratios2 = zeros(1,nPhases);

phasedVol = dissolvedVol.*lung_mask;
netPhase = sum(phasedVol(:));
for iPhase=1:nPhases
    phasedVec = netPhase*exp(1i*phaseVec(iPhase));
    ratios(iPhase) = real(phasedVec)/imag(phasedVec);
end

figure();
plot(phaseVec,ratios,'-b','LineWidth',3);


