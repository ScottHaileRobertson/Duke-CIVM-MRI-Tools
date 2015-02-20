function data = ScaleDataZeroToOne(data)
data = data-min(data(:));
data = data/max(data(:));