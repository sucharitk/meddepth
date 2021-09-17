function medriv_med_topos(plot_type, group, pthresh, nfreqs, grayscale,...
    do_fdr)


load('chanlocs.mat')
chanlocs = chanlocs32;
labels = {chanlocs.labels};

chanrem = {'TP9', 'TP10', 'FT10', 'FT9', 'FP1', 'Fp2'};
chanrem = get_channels_from_labels(labels, chanrem);

scale_mark = 1;
cbar = false;
labelposoffset = 1;

switch plot_type
    case 0
        if grayscale
            maplims = [-8 -4; 4 8; 0 5; 0 5];
        else
            %             maplims = [-8 0; 0 8; 0 5; 0 5];
            maplims = 8*[-1 0; 0 1; 0 5; 0 5];
        end
        fname = ['medpheno_alpha_theta_interaction_' group '.csv'];
        tplot_int(fname, pthresh, do_fdr, nfreqs, chanlocs, chanrem, ...
            scale_mark, cbar, grayscale,...
            labelposoffset, maplims)

    case 1
        if grayscale
            maplims = [-5 -2; 2 5; 2 5; 2 5];
        else
            maplims = [-1 0; 0 1; -1 1; -1 1]*6;
        end
        
        fname = ['medpheno_alpha_theta_interaction_' group '.csv'];
        tplot_int(fname, pthresh, do_fdr, nfreqs, chanlocs, chanrem, ...
            scale_mark, cbar, grayscale,...
            labelposoffset, maplims)
        
    case 2
        if grayscale
            maplims = [-5 -2; 2 5; 2 5; 2 5];
        else
            maplims = [-1 0; 0 1; -1 1; -1 1]*6;
        end
        postfix = '_last6';
        fname = ['medpheno_alpha_theta_interaction_' group postfix '.csv'];
        tplot_int(fname, pthresh, do_fdr, nfreqs, chanlocs, chanrem, ...
            scale_mark, cbar, grayscale,...
            labelposoffset, maplims)
        
    case 3
        if grayscale
            maplims = [-5 -2; 2 5; 2 5; 2 5];
        else
            maplims = [-1 0; 0 1; -1 1; -1 1]*6;
        end
        postfix = '_nochant';
        fname = ['medpheno_alpha_theta_interaction_' group postfix '.csv'];
        tplot_int(fname, pthresh, do_fdr, nfreqs, chanlocs, chanrem, ...
            scale_mark, cbar, grayscale,...
            labelposoffset, maplims)
        
    case 4
        if grayscale
            maplims = [-5 -2; 2 5; 2 5; 2 5];
        else
            maplims = [-1 0; 0 1; -1 1; -1 1]*6;
        end
        postfix = '_chantonly';
        fname = ['medpheno_alpha_theta_interaction_' group postfix '.csv'];
        tplot_int(fname, pthresh, do_fdr, nfreqs, chanlocs, chanrem, ...
            scale_mark, cbar, grayscale,...
            labelposoffset, maplims)
        
    case 5
        if grayscale
            maplims = [2 5];
        else
            maplims = [-1 1]*.015;
        end
        
        fname = ['medpheno_theta_' group '.csv'];
        
        cr = csvread(['medpheno_data/' fname]);
        cr2 = reshape(cr, [32 2]);
        cr2 = cr2(~chanrem, :);

        chanlocrem = chanlocs(~chanrem);
        
        frpow = (cr2(:,1));
        frpval = (cr2(:,2));
        
        if do_fdr
            frpval = fdr(frpval);
        end
        sigchans = find(frpval<pthresh);
        figure
        topoplot_sigp((frpow),chanlocrem,sigchans,maplims,'t',...
            scale_mark,labelposoffset,cbar,grayscale)
        
        
    case 6
        if grayscale
            maplims = [-5 5];
        else
            maplims = [-1 1]*.015;
        end
        
        fname = ['medpheno_gamma_' group '.csv'];
        
        cr = csvread(['medpheno_data/' fname]);
        cr2 = reshape(cr, [32 2]);
        cr2 = cr2(~chanrem, :);

        chanlocrem = chanlocs(~chanrem);
        
        frpow = (cr2(:,1));
        frpval = (cr2(:,2));
        
        if do_fdr
            frpval = fdr(frpval);
        end
        sigchans = find(frpval<pthresh);
        figure
        topoplot_sigp((frpow),chanlocrem,sigchans,maplims,'t',...
            scale_mark,labelposoffset,cbar,grayscale)
        
    case 7
        if grayscale
            maplims = [-5 5];
        else
            maplims = [-1 1]*.01;
        end
        
        fname = ['medpheno_beta_' group '.csv'];
        
        cr = csvread(['medpheno_data/' fname]);
        cr2 = reshape(cr, [32 2]);
        cr2 = cr2(~chanrem, :);

        chanlocrem = chanlocs(~chanrem);
        
        frpow = (cr2(:,1));
        frpval = (cr2(:,2));
        
        if do_fdr
            frpval = fdr(frpval);
        end
        sigchans = find(frpval<pthresh);
        figure
        topoplot_sigp((frpow),chanlocrem,sigchans,maplims,'t',...
            scale_mark,labelposoffset,cbar,grayscale)
        
    case 8
        if grayscale
            maplims = [-5 5];
        else
            maplims = [-2 2]*.01;
        end
        
        fname = ['medpheno_alpha_' group '.csv'];
        
        cr = csvread(['medpheno_data/' fname]);
        cr2 = reshape(cr, [32 2]);
        cr2 = cr2(~chanrem, :);

        chanlocrem = chanlocs(~chanrem);
        
        frpow = (cr2(:,1));
        frpval = (cr2(:,2));
        
        if do_fdr
            frpval = fdr(frpval);
        end
        sigchans = find(frpval<pthresh);
        figure
        topoplot_sigp((frpow),chanlocrem,sigchans,maplims,'t',...
            scale_mark,labelposoffset,cbar,grayscale)
        
    case 9
        if grayscale
            maplims = [-5 -2; 2 5; 2 5; 2 5];
        else
            maplims = [-1 0; 0 1; -1 1; -1 1]*6;
        end
        
        fname = ['medpheno_gamma_interaction_' group '.csv'];
        tplot_int(fname, pthresh, do_fdr, nfreqs, chanlocs, chanrem, ...
            scale_mark, cbar, grayscale,...
            labelposoffset, maplims)
        
    case 10
        if grayscale
            maplims = [-5 -2; 2 5; 2 5; 2 5];
        else
            maplims = [-1 0; 0 1; -1 1; -1 1]*5;
        end
        nfsign = [-1 1 -1 1];

        fname = 'medpheno_dl4_topo.csv';
        tplot_int(fname, pthresh, do_fdr, nfreqs, chanlocs, chanrem, ...
            scale_mark, cbar, grayscale,...
            labelposoffset, maplims, nfsign)

end
end

function tplot_int(fname, pthresh, do_fdr, nfreqs, chanlocs, ...
    chanrem, scale_mark, cbar, grayscale,...
            labelposoffset, maplims, nfsign)
     
cr = csvread(['medpheno_data/' fname]);
cr2 = reshape(cr, [32 nfreqs 2]);
cr2 = cr2(~chanrem, :, :);
nclust = 1;

chanlocrem = chanlocs(~chanrem);

% scalp_amp = sqrt(cr2(:,1,1));
% scalp2source(scalp_amp, chanlocrem, 2);

if ~exist('nfsign', 'var')
    nfsign = [-1 1 1 1];
end

figure
for nf = 1:nfreqs
    subplot(1,nfreqs,(nf-1)*nclust+1)
    frpow = nfsign(nf)*sqrt(cr2(:,nf,1));
    frpval = (cr2(:,nf,2));

    if do_fdr
        frpval = fdr(frpval);
    end
    sigchans = find(frpval<pthresh);

    
    topoplot_sigp((frpow),chanlocrem,sigchans,maplims(nf, :),'t',...
        scale_mark,labelposoffset,cbar,grayscale)
    
end
end

function tplot_main(fname, pthresh, do_fdr, nfreqs, chanlocs, chanrem, ...
    scale_mark, cbar, grayscale,...
            labelposoffset)
cr = csvread(fname);

cr2 = reshape(cr, [5 32 nfreqs 2]);
cr2 = cr2(:, ~chanrem, :, :);

nclust = 5;
maplims = [-4 4];

figure
for nf = 1:nfreqs
    for nc = 1:nclust
        subplot(nfreqs,5,(nf-1)*nclust+nc)
        frpow = (cr2(nc,:,nf,1));
        frpval = (cr2(nc,:,nf,2));
        %         frpval = (cr2(nc,:,nf,1));
        if do_fdr
            frpval = fdr(frpval);
        end
        sigchans = find(frpval<pthresh);
        %         sigchans = find(abs(frpval)>2.5);
        chanlocrem = chanlocs(~chanrem);

        topoplot_sigp((frpow),chanlocrem,sigchans,maplims,'t',...
            scale_mark,labelposoffset,cbar, grayscale)
        
    end
end
end