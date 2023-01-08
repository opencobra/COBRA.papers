curFolder = fileparts(mfilename('fullpath'));

%plotter = str2func('semilogx');
plotter = str2func('loglog');

%%%%%%%%%%%%%%%%%%%%%
%% 1. Rounding Time vs Dim
%subplot(1, 2, 1);
figure;
roundingTime(curFolder, plotter)

lgd = legend('BarrierRound', 'PolyRound', 'IsoRound', 'Location', 'northwest');
lgd.FontSize = 15;
ax = gca; 
ax.FontSize = 16; 

xlim([20 1e4]); %ylim([1e-5 1e5]);

%% 2. Mixing Time vs Dim
figure;
%subplot(1, 2, 2);
quantVSdim(curFolder, '2logBarrier/logBarrier_test', plotter, 'ro', "mixing")
quantVSdim(curFolder, '1polyRound/PolyRound_test', plotter, 'bo', "mixing")
quantVSdim(curFolder, '3isotropic/isotropic_test', plotter, 'go', "mixing")

lgd = legend('BarrierRound', 'PolyRound', 'IsoRound', 'Location', 'northwest');
lgd.FontSize = 15;
ax = gca; 
ax.FontSize = 16; 
xlim([10 1e3]); %ylim([10 1e9]);

function roundingTime(curFolder, plotter)
    %%%%%%%%%%%%
    % 2. Log-Barrier
    pathdir = 'Instances/2testbed-logBarrier';
    matfiles = dir(fullfile(fullfile(curFolder, pathdir), '*.mat'));

    numPoly = length(matfiles);
    dim = zeros(numPoly, 1);
    rtime = zeros(numPoly, 1);
    for idx = 1:numPoly
        load(strcat(matfiles(1).folder,'/',matfiles(idx).name));
        dim(idx) = P.dim;
        rtime(idx) = P.roundTime;
    end

    h = plotter(dim, rtime, 'ro');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    title('Rounding Time', 'FontSize', 18);
    xlabel('Dimension', 'FontSize', 18); ylabel('Time (sec)', 'FontSize', 18);
    hold on;
    fit = polyfit(log(dim), log(rtime), 1);
    newdim = dim;
    newdim = [20, 1100, 1000, 10000];
    z = polyval(fit, log(newdim));
    plotter(newdim, exp(z), 'r');
    algo = pathdir(1:4);
    fprintf(strcat(algo, ": rTime/Dim = %f & Coeff = %f\n"), fit(1), exp(fit(2)));

    % 1. PolyRound
    pathdir = '1polyRound/PolyRound/PolyRound/output';
    matfiles = csvread(fullfile(fullfile(curFolder, pathdir), 'roundTimes.csv'));
    rtime = matfiles(2:end,2);
    rtime = rtime(rtime>0);
    % only for p1~12, p15~18
    dim = [24, 58, 106, 191, 327, 430, 527, 580, 595, 607, 634, 932, 142, 544, 1056, 1410];

    h = plotter(dim, rtime, 'bo');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    title('Rounding Time', 'FontSize', 18);
    xlabel('Dimension', 'FontSize', 18); ylabel('Time (sec)', 'FontSize', 18);
    hold on;
    fit = polyfit(log(dim), log(rtime), 1);
    newdim = dim;
    newdim = [20, 1100, 1000, 10000];
    z = polyval(fit, log(newdim));
    plotter(newdim, exp(z), 'b');
    algo = pathdir(1:4);
    fprintf(strcat(algo, ": rTime/Dim = %f & Coeff = %f\n"), fit(1), exp(fit(2)));

    %%%%%%%%%%%%
    % 3. Isotropic
    pathdir = 'Instances/3testbed-isotropic';
    matfiles = dir(fullfile(fullfile(curFolder, pathdir), '*.mat'));

    numPoly = length(matfiles);
    dim = zeros(numPoly, 1);
    rtime = zeros(numPoly, 1);
    for idx = 1:numPoly
        load(strcat(matfiles(1).folder,'/',matfiles(idx).name));
        dim(idx) = size(P.A, 2);
        rtime(idx) = P.roundTime;
    end

    h = plotter(dim, rtime, 'go');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    title('Rounding Time', 'FontSize', 18);
    xlabel('Dimension', 'FontSize', 18); ylabel('Time (sec)', 'FontSize', 18);
    hold on;
    fit = polyfit(log(dim), log(rtime), 1);
    newdim = dim;
    newdim = [20, 1100, 1000, 10000];
    z = polyval(fit, log(newdim));
    plotter(newdim, exp(z), 'g');
    algo = pathdir(1:4);
    fprintf(strcat(algo, ": rTime/Dim = %f & Coeff = %f\n"), fit(1), exp(fit(2)));
end

function quantVSnnz(curFolder, pathdir, plotter, color, quant)
    matfiles = dir(fullfile(fullfile(curFolder, pathdir), '*.mat'));
    
    datapath = fullfile(curFolder, '/Instances/1chrr/', '*.mat');
    matfiles2 = dir(datapath);

    numPoly = length(matfiles);
    vnnz = zeros(numPoly, 1);
    time = zeros(numPoly, 1);
    step = zeros(numPoly, 1);
    for idx = 1:numPoly
        result = load(strcat(matfiles(1).folder,'/',matfiles(idx).name));
        if idx < 22
            load(fullfile(curFolder, 'Instances/1chrr/', matfiles2(idx).name));

            if result.exps.ess >= 10
                vnnz(idx) = nnz(P.A_eq);
                time(idx) = result.exps.sampleTime/result.exps.ess;
                step(idx) = result.exps.step/result.exps.ess;
            else
                fprintf("%s: Ess %d is too small\n", matfiles(idx).name, result.exps.ess) 
            end
        else
            if result.exps.ess >= 10
                vnnz(idx) = 295946;
                time(idx) = result.exps.sampleTime/result.exps.ess;
                step(idx) = result.exps.step/result.exps.ess;
            else
                fprintf("%s: Ess %d is too small\n", matfiles(idx).name, result.exps.ess) 
            end
        end
    end

    [vnnz, seq] = sort(vnnz); time = time(seq); step = step(seq);
    filter = vnnz>0;
    vnnz= vnnz(filter); time = time(filter); step = step(filter);
    
    if quant == "sample"
        h = plotter(vnnz, time, color);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        title('Sampling Time', 'FontSize', 15);
        xlabel('NNZ', 'FontSize', 15); ylabel('Time/Sample (s)', 'FontSize', 15);
        hold on;
        fit = polyfit(log(vnnz), log(time), 1);
        vnnz = [100, 1000, 10000, 100000, 1e6];
        z = polyval(fit, log(vnnz));
        plotter(vnnz, exp(z), color(1));
        
        algo = pathdir(1:4);
        fprintf(strcat(algo, ": Time/NNZ = %f\n"), fit(1));
    elseif quant == "mixing"
        h = plotter(vnnz, step, color);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        title('Mixing Rate', 'FontSize', 15);
        xlabel('NNZ', 'FontSize', 15); ylabel('Step/Sample', 'FontSize', 15);
        hold on;
        fit = polyfit(log(vnnz), log(step), 1);
        vnnz = [100, 1000, 10000, 100000, 1e6];
        z = polyval(fit, log(vnnz));
        plotter(vnnz, exp(z), color(1));
        
        algo = pathdir(1:4);
        fprintf(strcat(algo, ": Step/NNZ= %f\n"), fit(1));
    end
end

function quantVSdim(curFolder, pathdir, plotter, color, quant)
    matfiles = dir(fullfile(fullfile(curFolder, pathdir), '*.mat'));

    numPoly = length(matfiles);
    dim = zeros(numPoly, 1);
    time = zeros(numPoly, 1);
    step = zeros(numPoly, 1);
    for idx = 1:numPoly
        result = load(strcat(matfiles(1).folder,'/',matfiles(idx).name));

        if result.exps.ess >= 10
            dim(idx) = result.exps.dim;
            time(idx) = result.exps.sampleTime/result.exps.ess;
            step(idx) = result.exps.step/result.exps.ess;
        else
            fprintf("%s: Ess %d is too small\n", matfiles(idx).name, result.exps.ess) 
        end
    end

    [dim, seq] = sort(dim); time = time(seq); step = step(seq);
    filter = dim>0;
    dim = dim(filter); time = time(filter); step = step(filter);
    
    if quant == "sample"
        h = plotter(dim, time, color);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        title('Sampling Time', 'FontSize', 18);
        xlabel('Dimension', 'FontSize', 18); ylabel('Time/Sample (s)', 'FontSize', 18);
        hold on;
        fit = polyfit(log(dim), log(time), 1);
        newdim = dim;
        newdim = [10, 1000];
        z = polyval(fit, log(newdim));
        plotter(newdim, exp(z), color(1));
        algo = pathdir(1:4);
        fprintf(strcat(algo, ": Time/Dim = %f & Coeff = %f\n"), fit(1), exp(fit(2)));
    elseif quant == "mixing"
        h = plotter(dim, step, color);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        title('Mixing Rate', 'FontSize', 15);
        xlabel('Dimension', 'FontSize', 15); ylabel('Step/ESS', 'FontSize', 15);
        hold on;
        fit = polyfit(log(dim), log(step), 1);
        newdim = dim;
        newdim = [10, 1000];
        z = polyval(fit, log(newdim));
        plotter(newdim, exp(z), color(1));
        algo = pathdir(1:4);
        fprintf(strcat(algo, ": Step/Dim = %f & Coeff = %f\n"), fit(1), exp(fit(2)));
    end
end

function roundResult2(curFolder, pathdir, plotter, polytopes, color, dimOpt)
    matfiles = dir(fullfile(fullfile(curFolder, pathdir), '*.mat'));
    datapath = fullfile(curFolder, '/Instances/1chrr/', '*.mat');
    matfiles2 = dir(datapath);
    
    numPoly = length(matfiles2);
    dim = zeros(numPoly, 1);
    time = zeros(numPoly, 1);
    
    for idx = 1:numPoly
        load(fullfile(curFolder, 'Instances/1chrr/', matfiles2(idx).name));
        name = matfiles(idx).name(17:end-4); 
        dim(idx) = polytopes.(name).(dimOpt)(2);
        time(idx) = P.roundTime;
    end

    [dim, seq] = sort(dim); time = time(seq);
    filter = dim>0;
    dim = dim(filter); time = time(filter);
    h = plotter(dim, time, color);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    title('Preparation Time', 'FontSize', 15);
    xlabel('Dimension', 'FontSize', 15); ylabel('Time', 'FontSize', 15);
    hold on;
    fit = polyfit(log(dim), log(time), 1);
    newdim = [60, 100, 1000, 10000, 100000, 1e6];
    z = polyval(fit, log(newdim));
    plotter(newdim, exp(z), color(1));
    algo = pathdir(1:4);
    fprintf(strcat(algo, ": RoundTime/Dim = %f & Coeff = %f\n"), fit(1), exp(fit(2)));
end

function roundResult1(curFolder, pathdir, plotter, polytopes, color, dimOpt)
    matfiles = dir(fullfile(fullfile(curFolder, pathdir), '*.mat'));

    numPoly = length(matfiles);
    dim = zeros(numPoly, 1);
    time = zeros(numPoly, 1);
    for idx = 1:numPoly
        result = load(strcat(matfiles(1).folder,'/',matfiles(idx).name));
        name = matfiles(idx).name(17:end-4); 

        if result.exps.ess >= 10
            dim(idx) = polytopes.(name).(dimOpt)(2);
            time(idx) = result.exps.roundTime;
        else
            fprintf("%s: Ess %d is too small\n", matfiles(idx).name, result.exps.ess) 
        end
    end

    [dim, seq] = sort(dim); time = time(seq);
    filter = dim>0;
    dim = dim(filter); time = time(filter);
    h = plotter(dim, time, color);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    title('Preparation Time', 'FontSize', 15);
    xlabel('Dimension', 'FontSize', 15); ylabel('Time', 'FontSize', 15);
    hold on;
    fit = polyfit(log(dim), log(time), 1);
    newdim = [60, 100, 1000, 10000, 100000, 1e6];
    z = polyval(fit, log(newdim));
    plotter(newdim, exp(z), color(1));
    algo = pathdir(1:4);
    fprintf(strcat(algo, ": RoundTime/Dim = %f & Coeff = %f\n"), fit(1), exp(fit(2)));
end

function plotResult2(curFolder, pathdir, plotter, color)
    matfiles = dir(fullfile(fullfile(curFolder, pathdir), '*.mat'));
    
    datapath = fullfile(curFolder, '/Instances/1chrr/', '*.mat');
    matfiles2 = dir(datapath);

    numPoly = length(matfiles);
    vnnz = zeros(numPoly, 1);
    time = zeros(numPoly, 1);
    step = zeros(numPoly, 1);
    for idx = 1:numPoly
        result = load(strcat(matfiles(1).folder,'/',matfiles(idx).name));
        if idx < 22
            load(fullfile(curFolder, 'Instances/1chrr/', matfiles2(idx).name));

            if result.exps.ess >= 10
                vnnz(idx) = nnz(P.A_eq);
                time(idx) = result.exps.sampleTime/result.exps.ess;
                step(idx) = result.exps.step/result.exps.ess;
            else
                fprintf("%s: Ess %d is too small\n", matfiles(idx).name, result.exps.ess) 
            end
        else
            if result.exps.ess >= 10
                vnnz(idx) = 295946;
                time(idx) = result.exps.sampleTime/result.exps.ess;
                step(idx) = result.exps.step/result.exps.ess;
            else
                fprintf("%s: Ess %d is too small\n", matfiles(idx).name, result.exps.ess) 
            end
        end
    end

    [vnnz, seq] = sort(vnnz); time = time(seq); step = step(seq);
    filter = vnnz>0;
    vnnz= vnnz(filter); time = time(filter); step = step(filter);
    subplot(1,2,1);
    plotter(vnnz, time, color)
    title('Sampling Time', 'FontSize', 15);
    xlabel('NNZ', 'FontSize', 15); ylabel('Time/Sample (s)', 'FontSize', 15);
    hold on;

    subplot(1,2,2)
    plotter(vnnz, step, color)
    title('Mixing Time', 'FontSize', 15);
    xlabel('NNZ', 'FontSize', 15); ylabel('Step/Sample', 'FontSize', 15);
    hold on;

    fit1 = polyfit(log(vnnz), log(time), 1);
    fit2 = polyfit(log(vnnz), log(step), 1);
    
    algo = pathdir(1:4);
    fprintf(strcat(algo, ": Time/NNZ = %f, Step/NNZ= %f\n"), fit1(1), fit2(1));
end

function plotResult_RHMC(curFolder, pathdir, plotter, polytopes, color)
    matfiles = dir(fullfile(fullfile(curFolder, pathdir), '*.mat'));

    numPoly = length(matfiles);
    dim = zeros(numPoly, 1);
    time = zeros(numPoly, 1);
    step = zeros(numPoly, 1);
    for idx = 1:numPoly
        result = load(strcat(matfiles(1).folder,'/',matfiles(idx).name));
        name = matfiles(idx).name(17:end-4); 

        if result.exps.ess >= 10
            dim(idx) = polytopes.(name).processedSize(2);
            time(idx) = result.exps.sampleTime/result.exps.ess;
            step(idx) = result.exps.step/result.exps.ess;
        else
            fprintf("%s: Ess %d is too small\n", matfiles(idx).name, result.exps.ess) 
        end
    end

    [dim, seq] = sort(dim); time = time(seq); step = step(seq);
    subplot(1,2,1);
    plotter(dim, time, color)
    title('Sampling Time', 'FontSize', 15);
    xlabel('Dimension', 'FontSize', 15); ylabel('Time/Sample (s)', 'FontSize', 15);
    hold on;

    subplot(1,2,2)
    plotter(dim, step, color)
    title('Mixing Time', 'FontSize', 15);
    xlabel('Dimension', 'FontSize', 15); ylabel('Step/Sample', 'FontSize', 15);
    hold on;

    fit1 = polyfit(log(dim), log(time), 1);
    fit2 = polyfit(log(dim), log(step), 1);
    
    algo = pathdir(1:4);
    fprintf(strcat(algo, ": Time/Dim = %f, Step/Dim = %f\n"), fit1(1), fit2(1));
end

function plotResult(curFolder, pathdir, plotter, polytopes, color, dimOpt)
    matfiles = dir(fullfile(fullfile(curFolder, pathdir), '*.mat'));

    numPoly = length(matfiles);
    dim = zeros(numPoly, 1);
    time = zeros(numPoly, 1);
    step = zeros(numPoly, 1);
    for idx = 1:numPoly
        result = load(strcat(matfiles(1).folder,'/',matfiles(idx).name));
        name = matfiles(idx).name(17:end-4); 

        if result.exps.ess >= 10
            dim(idx) = polytopes.(name).(dimOpt)(2);
            time(idx) = result.exps.sampleTime/result.exps.ess;
            step(idx) = result.exps.step/result.exps.ess;
        else
            fprintf("%s: Ess %d is too small\n", matfiles(idx).name, result.exps.ess) 
        end
    end

    [dim, seq] = sort(dim); time = time(seq); step = step(seq);
    filter = dim>0;
    dim = dim(filter); time = time(filter); step = step(filter);
    subplot(1,2,1);
    plotter(dim, time, color)
    title('Sampling Time', 'FontSize', 15);
    xlabel('Dimension', 'FontSize', 15); ylabel('Time/Sample (s)', 'FontSize', 15);
    hold on;

    subplot(1,2,2)
    plotter(dim, step, color)
    title('Mixing Time', 'FontSize', 15);
    xlabel('Dimension', 'FontSize', 15); ylabel('Step/Sample', 'FontSize', 15);
    hold on;

    fit1 = polyfit(log(dim), log(time), 1);
    fit2 = polyfit(log(dim), log(step), 1);
    
    algo = pathdir(1:4);
    fprintf(strcat(algo, ": Time/Dim = %f, Step/Dim = %f\n"), fit1(1), fit2(1));
end

