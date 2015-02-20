function [ structToOverride ] = overrideStructFields( structToOverride, overridingStruct)
%OVERRIDESTRUCTFIELDS Overrides fields of a structure with fields of
%another struct

% Get list of all field names in overriding struct
fieldVals = fieldnames(overridingStruct)
nField = length(fieldVals);

% Loop through fields and update them
for iField = 1:nField
	newStructVal = getfield(overridingStruct,fieldVals{iField});
	if(isfield(structToOverride, fieldVals{iField}))
		oldStructVal = getfield(structToOverride,fieldVals{iField});
		
		% If there are substructures, we need to recurse
		if(isstruct(oldStructVal))
			if(isstruct(newStructVal))
				% Recurse
				newStruct = overrideStructFields(oldStructVal,newStructVal);
				structToOverride = setfield(structToOverride,fieldVals{iField},newStruct);
			else
				error('Both values must be structures');
			end
		else
			% Otherwise we are at a node and can set the field value
			if(~isstruct(newStructVal))
				structToOverride = setfield(structToOverride,fieldVals{iField},newStructVal);
			else
				error('Both values must not be structures');
			end
		end
	else
		warning(['Adding field ' fieldVals{iField} ' to struct']);
		structToOverride = setfield(structToOverride,fieldVals{iField},newStructVal);
	end
end

end

