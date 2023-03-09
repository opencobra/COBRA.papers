%problemlist = {'LPnetlib@"lp_25fv47"'};%, 'LPnetlib@"lp_80bau3b"', 'LPnetlib@"lp_afiro"', 'LPnetlib@"lp_agg"', 'LPnetlib@"lp_beaconfd"', 'LPnetlib@"lp_blend"', 'LPnetlib@"lp_degen2"', 'LPnetlib@"lp_degen3"', 'LPnetlib@"lp_etamacro"', 'LPnetlib@"lp_scorpion"', 'LPnetlib@"lp_sierra"', 'LPnetlib@"lp_truss"'};
%problemlist = {'metabolic/Recon1'};
%problemlist = {'tvball/tvball'};
result = struct;
attempts = 10000; skip = 100;
john_ellip = 0; roundparam = 1; 
refSamples = 500000; % set it if roundparam=2

trials = 5;
tlist = zeros(trials, 1); slist = zeros(trials, 1);
%problem = problemTVballDifficult(100);
problem = loadProblem_rhmc('test');
[P, ~, ~] = fullify_and_sample(problem, 1, 1, john_ellip, roundparam, refSamples);

fprintf('Parameters - John: %d, Round: %d, Repeat: %d\n', john_ellip, roundparam, trials);
for i = 1:trials
    [P, o, sampleTime] = fullify_and_sample(problem, attempts, skip, john_ellip, roundparam, refSamples, P); 
    %disp(size(P.A));

    nSamples = size(o, 2);
    tlist(i) = sampleTime / nSamples;
    slist(i) = attempts*skip / nSamples;
    
    %test.thinned = true; test.samples = o; test.dim = size(problem.Aeq, 2) - size(problem.Aeq, 1);
    %pVal = uniformTest(test, problem, struct('toPlot', true));
    
    %% Check feasibility of Samples
    assert(min(o-problem.lb, [], 'all') >= 0)
    assert(min(problem.ub-o, [], 'all') >= 0)
    assert(norm(problem.beq - problem.Aeq * o, Inf) <= 1e-8)
end
fprintf('Eff Sampling Time : %e\n', mean(tlist))
fprintf('# Eff Steps: %e\n', mean(slist))

% for i = 1:length(problemlist)
%     %problem = loadProblem_rhmc(problemlist{i});
%     problem = problemTVballDifficult(100);
%     %problem = loadProblem(problemlist{i});
%     [P, o, sampleTime] = fullify_and_Sample(problem, 0, attempts, skip);
%     disp(size(P.A));
%     
%     %[o, sampleTime] = Sample_time(attempts, skip, P, [], 0.1);
%     nSamples = size(o, 2);
%     stimes = sampleTime / nSamples;
%     steps = attempts*skip / nSamples;
%     
%     exps = struct; exps.size = size(P.A); 
%     exps.stimes = stimes; exps.steps = steps;
%     
%     tmp = split(problemlist{i}, '/');
%     fieldname = tmp{2};
%     result.(fieldname) = exps;
%     disp([problemlist{i} ' done']);
% end
% 
% for i = 1:length(problemlist)
%     tmp = split(problemlist{i}, '/');
%     fieldname = tmp{2};
%     disp(fieldname);
%     disp(result.(fieldname));
% end