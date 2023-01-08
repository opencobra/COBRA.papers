function [P, o, sampleTime] = fullify_and_sample(problem, attempts, skips, john_ellip, roundparam, refSamples, P)
      if nargin < 7
          %% prepare input polytope P={x: Ax=b, l<x<u}
          % Run John's ellipspod in Reflection walk (depending on the parameter
          % 'john_elip')
          refwalk = sample(problem, 1, struct('presolve', john_ellip));
          P.x0 = refwalk.prepare.x0; dim = length(P.x0);
          P.idx = refwalk.prepare.idx; P.scale = refwalk.prepare.scale;

          % rounding based on samples from Reflection walk
          refwalk = sample(problem, refSamples, struct('presolve', false, 'prepare', refwalk.prepare));
          [P.R, P.R_shift] = round_refwalk(refwalk.samples);
          P.center = zeros(size(P.R, 2), 1);
          P.A_eq = refwalk.prepare.A(:, 1:dim) * P.R;
          P.b_eq = refwalk.prepare.b - refwalk.prepare.A(:, 1:dim) * P.R_shift;
          P.A = [P.R; -P.R];
          P.b = [double(refwalk.prepare.ub(1:dim)-P.R_shift); -double(refwalk.prepare.lb(1:dim)-P.R_shift)];
      end
      
      %% CHAR sampling
      tic;
      o = Sample_time(attempts, skips, P, [], 0.1);
      
      %% Recover samples 
      % transform back 1 (to go back to Ax=b space)
      if roundparam == 2
          o = P.R*o + P.R_shift;
      else
          o = P.N*o + P.p_shift;
      end
      
      % transform back 2 (due to preprocessing at the beginning)
      if john_ellip
          samplesNew = P.x0 * ones(1, size(o, 2));
          samplesNew(P.idx,:) = samplesNew(P.idx,:) + P.scale.*o;
          o = samplesNew;
      end
      
      % consider only effective samples
      ess = effectiveSampleSize(o);
      nSamples = min(ess);
      gap = round(size(o,2) / nSamples);
      o = o(:,gap:gap:end);
   
      sampleTime = toc;
end

function [T, T_shift] = round_refwalk(o)
    nSamples = size(o, 2);
    T_shift = mean(o, 2);
    o_origin = o - T_shift;
    [U, S, ~] = svd(o_origin);
    s = diag(S); s = s(s>1e-10); degenDim = length(s);
    S = diag(s); U = U(:, 1:degenDim);
    T = U*S/sqrt(nSamples-1);
end