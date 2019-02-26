function obj = wf_catb(shape_name,Fy,design_basis)

if nargin < 3
    design_basis = 'AISC2016';
end

switch design_basis
    case 'AISC2010'
        obj = wf_catb_AISC2010(shape_name,Fy);
    case 'AISC2016'
        obj = wf_catb_AISC2016(shape_name,Fy);
    otherwise
        error('Unknown design basis: %s',design_basis)
end

end
