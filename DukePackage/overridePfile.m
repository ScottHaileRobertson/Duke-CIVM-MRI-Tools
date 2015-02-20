function pfile = overridePfileHeader(pfile, pfileOverride)
% Override rdb
if(~isempty(pfileOverride.rdb))
    pfile.rdb = overrideStructFields(pfile.rdb, pfileOverride.rdb);
end

% Override exam
if(~isempty(pfileOverride.exam))
    pfile.exam = overrideStructFields(pfile.exam, pfileOverride.exam);
end

% Override series
if(~isempty(pfileOverride.series))
    pfile.series = overrideStructFields(pfile.series, pfileOverride.series);
end

% Override image
if(~isempty(pfileOverride.image))
    pfile.image = overrideStructFields(pfile.image, pfileOverride.iamge);
end
end