%VOLUME This object just holds an N-dimmensional volume (a N-D space) and
%       some basic properties of the volume. The benefit to this class is
%       that its an object, so it can be passed by refference to functions
%       (saves memory).
%
%   Copyright: 2012 Scott Haile Robertson.
%   Website: www.ScottHaileRobertson.com
%   $Revision: 1.1 $  $Date: 2012/04/04 $
%       4/5/2013 Added dimmension names
classdef LinkedVolume %< Volume
    properties
        Volumes;
    end
    
    methods
        function this = LinkedVolume(varargin)
            if(nargin > 1)
                if(iscell(varargin{1}))
                    vols = varargin{1};
                    for i=1:length(vols)
                        if(~isa(vols(i),'Volume'))
                            error('Linked Volume requires all Volumes as an input.');
                        end
                    end
                    this.Volumes = vols;
                else
                    error('Volumes must be passed in as a cell.');
                end
            end
        end
        
        function this = updateData(this, newData, volIdx)
            this.updateData(newData,1);
        end
        
        function this = updateData(this, newData, volIdx)
            this.Volumes.updateData(newData);
        end
        
        function this = setDimmensionNames(this, newDimNames)
            if(iscell(newDimNames))
                if(length(newDimNames)==this.NDims)
                    for i=1:this.NDims
                        if(~isa(newDimNames{i},'char'))
                            if(isnumeric(newDimNames{i}))
                                newDimNames{i} = num2str(newDimNames{i});
                            else
                                error('Dimmension names must be strings or numbers.');
                            end
                        end
                    end
                    this.DimNames = newDimNames;
                else
                    error(['Number of dimmension names must equal number of '...
                        'dimmensions in volume']);
                end
            else
                error('Dimmension names must be passed in as a cell.');
            end
        end
    end
end