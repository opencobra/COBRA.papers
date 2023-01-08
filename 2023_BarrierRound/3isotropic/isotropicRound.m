function o = isotropicRound(problem)
%Input: a structure P with the following fields
%  .Aeq
%  .beq
%  .lb
%  .ub
% describing the polytope {Aeq x = beq, lb <= x <= ub}
%Output:
% o - problem structure

% Run CRHMC in COBRA package to get O*(n) samples with transformation (that
% will be used to recover original coordinates)
opts = default_options();
opts.rawOutput = true;
opts.maxTime = +Inf;
roughDim = size(problem.Aeq, 2) - size(problem.Aeq, 1);
o = sample(problem, roughDim, opts);
degenDim = size(o.problem.A, 2) - size(o.problem.A, 1);

% Thin samples here to ensure independence
ess = min(effective_sample_size(o.chains), [], 1);
out = [];
N = size(o.chains, 3);
for i = 1:numel(ess)
    gap = ceil(N/ess(i));
    out_i = o.chains(i, :, 1:gap:N);
    out_i = reshape(out_i, [size(out_i,2) size(out_i,3)]);
    out = [out out_i];
end

% Compute isotropic transformation
[R, R_shift] = isoTransform(out, degenDim);
%nSamples = size(out, 2);
%[R, Rinv, R_shift] = isoTransform2(out(:, 1:round(nSamples/2)), nullDim);

% Check isotropy
%checkIsotropy(out(:, round(nSamples/2)+1:end), Rinv, R_shift);

% Full-dimensional description Ax<=b
A = [R; -R];
b = [o.problem.barrier.ub'-R_shift; -(o.problem.barrier.lb'-R_shift)];
center = zeros(size(A, 2), 1);

% Transformation for recovery
T = o.problem.T * R;
y = o.problem.T * R_shift + o.problem.y; 

o = struct;
o.A = A; o.b = b; o.x = center; % rounded polytope with feasible point x
o.T = T; o.y = y; % used to recover sample in original space by x' = Tx + y
end

function [T, T_shift] = isoTransform(o, nullDim)
    nSamples = size(o, 2);
    T_shift = mean(o, 2);
    o_origin = o - T_shift;
    [U, S, ~] = svd(o_origin);
    S = S(1:nullDim, 1:nullDim);
    U = U(:,1:nullDim);
    T = U*S/sqrt(nSamples);
end

function [T, Tinv, T_shift] = isoTransform2(o, nullDim)
    nSamples = size(o, 2);
    T_shift = mean(o, 2);
    o_origin = o - T_shift;
    [U, S, ~] = svd(o_origin);
    S = S(1:nullDim, 1:nullDim);
    U = U(:,1:nullDim);
    T = U*S/sqrt(nSamples);
    Tinv = sqrt(nSamples)*U'./diag(S);
end

function checkIsotropy(o, Rinv, R_shift)
    newSamples = Rinv*(o-R_shift);
    empCov = newSamples*newSamples'/size(o, 2);

    e = eig(empCov);
    plot(e, 'o')
end