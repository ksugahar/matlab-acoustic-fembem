function model = volFemBem(volFile, options)
%VOLFEMBEM Gypsilab-style entry point for .vol FEM/BEM prototypes.
%
% Example:
%   m  = volFemBem("mesh.vol");
%   uh = h1(m);       % H1 P1 tetrahedra
%   ah = hcurl(m);    % HCurl Nedelec0 tetrahedra
%   jh = rwg(m);      % boundary RWG triangle-pair dofs

arguments
    volFile (1,1) string
    options.GypsilabRoot (1,1) string = ""
end

model = volFemBemModel(volFile, options.GypsilabRoot);
model = addFirstOrderFemBemSpaces(model);
model.topology = buildFirstOrderTopology(model);

model.h1 = h1(model);
model.hcurl = hcurl(model);
model.rwg = rwg(model);
model.status = "vol_ready_first_order_h1_hcurl_rwg";
end
