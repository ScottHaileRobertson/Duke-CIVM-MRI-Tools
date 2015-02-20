function dc_sample_idx = calcDCsample(radialDistance)
if(any(radialDistance == 0))
	% Fit the last DC point
	dc_sample_idx = find(radialDistance==0,1,'last');
else
	% If no DC point, just fit the first point
	warning('No DC point in radialDistance: calculating flip angle from first sample...');
	dc_sample_idx = 1;
end