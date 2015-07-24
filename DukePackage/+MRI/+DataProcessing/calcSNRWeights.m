function 	weights = calcSNRWeights(pfile, dc_sample_idx, weight_type)

% Compensate for Decay via amplification
switch(weight_type)
	case 0
		% Apply uniform weights
		weights = ones(size(pfile.data));
	case 1
		% Apply frame by frame weighting based on DC signal
		weights = abs(repmat(pfile.data(dc_sample_idx,:),[size(pfile.data,1) 1]));
	otherwise
		error('Weight type not supported');
end
