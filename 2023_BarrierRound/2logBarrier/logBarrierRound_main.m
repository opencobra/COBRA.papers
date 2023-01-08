curFolder = fileparts(mfilename('fullpath'));
datapath = fullfile(fileparts(curFolder), '/Instances/0raw', '*.mat');
matfiles = dir(datapath);
numPoly = size(matfiles, 1);

for c = 1:numPoly
    model = load(fullfile(fileparts(curFolder), 'Instances/0raw', matfiles(c).name));
    model = cellfun(@(x)(model.(x)), fieldnames(model));
    P = struct;
    
    if isfield(model, 'S') % Bio
        P.Aeq = model.S;
        P.beq = model.b;
        P.lb = model.lb;
        P.ub = model.ub;
        if isfield(model,'dsense')
           I = (model.dsense == 'E');
           P.Aeq = [P.Aeq; model.C(I,:)];
           P.beq = [P.beq; model.d(I)];
           P.Aineq = model.C(~I,:);
           P.bineq = model.d(~I,:);
           flip = 1-2*(model.dsense(~I) == 'G');
           P.Aineq = flip.*P.Aineq;
           P.bineq = flip.*P.bineq;
        end
        P = standardize_problem(P);

        % Change to the form of Ax = b, l <= x <= u
        nP = size(P.Aeq, 2);
        nIneq = size(P.Aineq, 1);
        nEq = size(P.Aeq, 1);
         
        P.Aeq = [P.Aeq sparse(nEq, nIneq); P.Aineq speye(nIneq)];
        P.beq = [P.beq; P.bineq];
        P.lb = [P.lb; zeros(nIneq, 1)];
        P.ub = [P.ub; Inf*ones(nIneq, 1)];
    else % Net Lib
        P.Aeq = model.A;
        P.beq = model.b;
        P.lb = model.aux.lo;
        P.ub = model.aux.hi;
    end
    
    tic;
    P = logBarrierRound(P); % resulting in a constraint-based rounded polytope
    roundTime = toc;
    dim = size(P.P.A,2) - size(P.P.A,1);
    fprintf("Loading complete (%d/%d): %s\n", c, numPoly, matfiles(c).name)
    
    % Find a full-dimensional description of this degenerate subspace
    skip = 1;
    if size(P.P.A,2)-size(P.P.A,1)<1100
        P = fullification(P);
        skip = 0; 
    end
    P.dim = dim;
    P.roundTime = roundTime;
    P.skip = skip;
    save(fullfile(fileparts(curFolder), '/Instances/2testbed-logBarrier/', matfiles(c).name), 'P')
%     if size(P.A, 2) < 6000
%         save(fullfile(fileparts(curFolder), '/Instances/2testbed-logBarrier/', matfiles(c).name), 'P')
%     else
%         save(fullfile(fileparts(curFolder), '/Instances/2testbed-logBarrier/', matfiles(c).name), 'P', '-v7.3')
%     end
end