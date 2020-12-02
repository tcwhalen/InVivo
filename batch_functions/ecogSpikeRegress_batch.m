function [ data ] = ecogSpikeRegress_batch( data )
% Tim Whalen, last edited Nov 2020
% Batch processing for ecogSpikeRegress, see for details
% Denoises ECoG then performs ecog-spike regression on all files in data
% and plots regression coefficients and error
%
% data must include: ts, T, ecogs, ecog_noise_thresh, files_trunc (for title, only if plotting)
% To specify parameters, add the sub-struct "ecog_reg" to your data struct 
% containg:any of the parameters in ecogSpikeRegress, and/or the following:
%   min_rate: minimum rate for unit to be analyzed, in Hz (default = 5)
%   to_plot: 1 to plot regression coeffs and error. Autocorrelation is also
%       plotted iff data has osc struct as outputted by renewalPSD_phaseShift 
%       and osc.min_rate == ecog_reg.min_rate. Will plot all units (whether
%       or not they have osc) regardless
%
% All outputs are the same as ecogSpikeRegress, wrapped in (nfiles,1) cells
% or vectors

% then get parameters or set to defaults
names = {'bin', 'step', 'ar_lag', 'nlag_past', 'nlag_fut', 'alpha', 'min_rate', 'do_sdf', 'sdf_std', 'backwards', 'to_plot'};
defaults = {.01, .01, 10, 100, 0, 0.05, 0, 1, 0.1, 0, 0};
[data, bin, step, ar_lag, nlag_past, nlag_fut, alpha, min_rate, do_sdf, sdf_std, backwards, to_plot] = extractDataInputs(data,'ecog_reg',names,defaults);

nfiles = data.nfiles;
data.ecog_reg.sig_lags = cell(nfiles,1);
data.ecog_reg.bs = cell(nfiles,1);
data.ecog_reg.ssrs = cell(nfiles,1);
data.ecog_reg.rmse_scaled = cell(nfiles,1);
data.ecog_reg.pttype = cell(nfiles,1);
data.ecog_reg.mean_amp = zeros(nfiles,1);
for f = 1:data.nfiles
    disp(['Computing ECoG regression model of file ' int2str(f) ' of ' int2str(nfiles)])
    ecog = data.ecogs{f};
    ecog = ecogThreshold(ecog,data.ecog_noise_thresh(f));
    ts = data.ts{f}(data.rates{f}>min_rate);
    [ecog_reg] = ecogSpikeRegress(ts, ecog, data.T(f), data.ecog_start(f), data.ECOG_FS(f), 'bin', bin, 'step', step, 'ar_lag', ar_lag, 'nlag_past', nlag_past, 'nlag_fut', nlag_fut, 'alpha', alpha, 'do_sdf', do_sdf, 'sdf_std', sdf_std, 'backwards', backwards);

    data.ecog_reg.sig_lags{f} = ecog_reg.sig_lags;
    data.ecog_reg.bs{f} = ecog_reg.bs;
    data.ecog_reg.ssrs{f} = ecog_reg.ssrs;
    data.ecog_reg.rmse_scaled{f} = ecog_reg.rmse_scaled;
    data.ecog_reg.pttype{f} = ecog_reg.pttype;
    data.ecog_reg.mean_amp(f) = ecog_reg.mean_amp;
end

if to_plot
    disp('Beginning plots')
    plot_acorr = 1;
    if ~isfield(data.osc,'acorr')
        disp('No autocorrelation found in data.osc. Skipping acorr plots')
        plot_acorr = 0;
    elseif data.osc.min_rate ~= data.ecog_reg.min_rate
        disp('Different min_rate for osc and ecog_reg, indices will be mismatched. Skipping acorr plots')
        plot_acorr = 0;
    end
    for f = 1:nfiles
        for u = 1:size(data.ecog_reg.bs{f},1)
            figure
            if plot_acorr
                subplot(3,1,1)
                plot(data.osc.acorrtimes,data.osc.acorr{f}(u,1:length(data.osc.acorrtimes)),'k')
                xlabel('Time (msec)')
                ylabel('Autocorr.')
                title('Spike Autocorrelation')
            end
            subplot(3,1,2)
            plot(1000*(-nlag_past*step:step:nlag_fut*step),data.ecog_reg.bs{f}(u,:))
            title('Spike-ECoG Regression Coefficients')
            ylabel('Regression Coeffs.')
            xlabel('<-- Past Spikes                 Lag (msec)               Future Spikes -->')
            subplot(3,1,3)
            plot(1000*(-nlag_past*step:step:nlag_fut*step),data.ecog_reg.ssrs{f}(u,:))
            title('Spike-ECoG Regression Residuals')
            ylabel('\Sigma (residuals^2)')
            xlabel('<-- Past Spikes                 Lag (msec)               Future Spikes -->')
            lag = data.ecog_reg.sig_lags{f}(u);
            if ~isnan(lag)
                hold on
                plot(1000*lag, min(data.ecog_reg.ssrs{f}(u,:)), 'r.','MarkerSize',10)
                text(1000*(lag+.02),min(data.ecog_reg.ssrs{f}(u,:)),int2str(round(lag*1000)),'Color','r','FontSize',10)
            end
            suptitle(['SNr-ECoG Regression (' data.files_trunc(f,:) ' Unit ' int2str(u) ')'])
            set(gcf, 'Position', [100, 800, 600 750])
        end
    end
end
end