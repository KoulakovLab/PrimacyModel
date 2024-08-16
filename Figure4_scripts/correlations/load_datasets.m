function [corr_mx, corr_mz, DX, DZ] = load_datasets(varargin)

load('../processed_connectivity.mat', 'processed_connectivity');

corr_mx = processed_connectivity.FlyEM.corr;
corr_mz = processed_connectivity.FAFB.corr;

DX = processed_connectivity.FlyEM.dist;
DZ = processed_connectivity.FAFB.dist;

if nargin > 0
    if strcmp(varargin{1}, 'shuffle')
        corr_mx = processed_connectivity.FlyEM.null_corr;
        corr_mz = processed_connectivity.FAFB.null_corr;
        
        DX = processed_connectivity.FlyEM.null_dist;
        DZ = processed_connectivity.FAFB.null_dist;
    end
end

end
