# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.2'
#       jupytext_version: 1.1.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# ## Setup

# %%
# %run setup.ipynb
rcParams['savefig.dpi'] = 200

# %%
genome = phase2_ar1.genome_agamp3

# %%
region_vgsc = SeqFeature('2L', 2358158, 2431617, genome=genome)

# %%
len(region_vgsc)

# %%
phase2_ar1.load_geneset_agamp44(attributes=['ID', 'Parent'])
geneset_agamp44 = geneset_to_pandas(phase2_ar1.geneset_agamp44)
geneset_agamp44_vgsc = geneset_agamp44.query(region_vgsc.query).copy()
geneset_davies = geneset_to_pandas(allel.FeatureTable.from_gff3('davies_vgsc_model.gff3', attributes=['ID', 'Parent']))
geneset_vgsc_combined = pandas.concat([geneset_agamp44_vgsc, geneset_davies])

# %%
tbl_davies_exons = (
    etl
    .fromdataframe(geneset_davies)
    .eq('type', 'CDS')
    .cutout('Parent', 'source', 'type', 'score', 'strand', 'phase')
    .merge(key=('start', 'end'))
    .movefield('seqid', 0)
)
tbl_davies_exons.displayall()

# %%
samples = pandas.read_csv(phase2_ar1.samples_fn, sep='\t')
samples.head()

# %%
sample_ids = samples.ox_code.values.tolist()
len(sample_ids)


# %% [markdown]
# ## Plot exons

# %%
def plot_exons(features, start, end, ax=None, xytext=None, y=0, height=1, gene_ec='#aaaaaa', gene_fc='none', exon_color='k', label_position='top'):
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 1))
        
    # plot whole gene span
    xmin = features.start.min()
    xmax = features.end.max()
    patch = plt.Rectangle((xmin, y), xmax-xmin+1, height, lw=2, edgecolor=gene_ec, facecolor=gene_fc)
    ax.add_patch(patch)
    
    for _, feature in features.iterrows():

        if feature.end >= start and feature.start <= end:
            x = feature.start
            width = feature.end - feature.start + 1
            patch = plt.Rectangle((x, y), width, height, lw=.2, edgecolor='k', facecolor='k')
            ax.add_patch(patch)
            
            try:
                xyt = xytext[feature.ID]
            except:
                xyt = (0, 8)
            if label_position == 'top':
                xy = (x + width/2, y + height)
                va = 'bottom'
                # leave xyt
            else:
                xy = (x + width/2, y)
                xyt = xyt[0], -xyt[1]
                va = 'top'
            ax.annotate(feature.ID, xy=xy, xytext=xyt, textcoords='offset points', ha='center', va=va, 
                        arrowprops=dict(arrowstyle='-', shrinkA=0, shrinkB=0, lw=.5, color='grey'))

    ax.set_xlim(start, end)
    ax.set_ylim(0, 1)
    
    
def plot_davies_exons(start, end, y=0, height=1, ax=None, label_position='top'):
    df = tbl_davies_exons.select('ID', lambda v: '-' not in v).todataframe()
    plot_exons(df, start=start, end=end, ax=ax, label_position=label_position,
               xytext={'7': (-5, 8),
                       '9': (3, 8),
                       '10': (4, 8),
                       '14': (3, 8),
                       '15': (-5, 8),
                       '17': (5, 8),
                       '18b+': (-7, 8),
                       '20c': (5, 8),
                       '20d': (-5, 8),
                       '21': (-4, 8),
                       '23f+': (-4, 8),
                       '24h+': (0, 16),
                       '25': (2, 8),
                       '26': (7, 16),
                       '27k': (5, 8),
                       '28': (-6, 16),
                       '29': (-1, 8),
                       '30': (0, 16),
                       '31': (2, 8),
                       '32': (4, 16),
                        },
               y=y, height=height)
    ax.set_ylim(-.1, 1.1)
    ax.set_yticks([])
    labels = 'exon', 
    handles = [plt.Rectangle((0, 0), 1, 1, color='k')]
    ax.legend(handles=handles, labels=labels, loc='center right', bbox_to_anchor=(0, .5))



# %%
fig, ax = plt.subplots(figsize=(15, 1))
sns.despine(ax=ax, left=True)
plot_davies_exons(start=region_vgsc.start - 5000, end=region_vgsc.end + 5000, ax=ax)


# %%
fig, ax = plt.subplots(figsize=(15, 1))
sns.despine(ax=ax, left=True, bottom=True)
plot_davies_exons(start=region_vgsc.start - 5000, end=region_vgsc.end + 5000, ax=ax, label_position='bottom')
ax.set_xticks([]);

# %%
fig, ax = plt.subplots(figsize=(15, 1))
sns.despine(ax=ax, left=True)
plot_davies_exons(start=region_vgsc.start - 5000, end=region_vgsc.start + 5000, ax=ax)


# %% [markdown]
# ## Plot sequence composition

# %%
def plot_seq_composition(genome, seqid, start, end, n_bins=500, ax=None, colors=None):

    # obtain reference sequence as numpy char array
    seq = np.asarray(genome[seqid])[start-1:end]

    # convert to lower-case
    seq = np.char.lower(seq)

    # locate nucleotides
    is_a = seq == b'a'
    is_c = seq == b'c'
    is_g = seq == b'g'
    is_t = seq == b't'
    is_n = seq == b'n'
    is_other = ~is_a & ~is_c & ~is_g & ~is_t & ~is_n
    assert np.sum(is_other) == 0

    # construct bins
    bins = np.linspace(0, len(seq), n_bins).astype(int)

    # count nucleotides
    h_a, _ = np.histogram(np.nonzero(is_a)[0], bins=bins)
    h_c, _ = np.histogram(np.nonzero(is_c)[0], bins=bins)
    h_g, _ = np.histogram(np.nonzero(is_g)[0], bins=bins)
    h_t, _ = np.histogram(np.nonzero(is_t)[0], bins=bins)
    h_n, _ = np.histogram(np.nonzero(is_n)[0], bins=bins)

    # plot
    left = bins[:-1] + start
    bottom = 0
    width = np.diff(bins)
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 1))
    if colors is None:
        palette = sns.color_palette('colorblind')
        colors = [palette[i] for i in [1, 0, 2, 4]] + ['k']
    else:
        assert len(colors) == 5, 'bad colors'
    for h, c, l in zip([h_a, h_t, h_g, h_c, h_n], colors, 'ATGCN'):
        ax.bar(left, h, width=width, bottom=bottom, color=c, align='edge', label=l)
        bottom += h
    ax.set_xlim(start, end)
    ax.set_yticks(ax.get_ylim())
    ax.set_yticklabels(['0%', '100%'])
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in ax.get_xticks()])
    
    # legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], loc='center left', bbox_to_anchor=(1, .5), prop=dict(family='monospace'), borderpad=0, ncol=2)
    


# %% [markdown]
# ## Plot accessibility

# %%
sns.palplot(sns.color_palette('Set3'))


# %%
def plot_accessibility(accessibility, seqid, start, end, n_bins, ax=None, colors=None):
 
    # extract accessibility map
    is_accessible = accessibility[seqid]['is_accessible'][start-1:end]

    # construct bins
    bins = np.linspace(0, len(is_accessible), n_bins).astype(int)

    # count bases
    h_a, _ = np.histogram(np.nonzero(is_accessible)[0], bins=bins)
    h_n, _ = np.histogram(np.nonzero(~is_accessible)[0], bins=bins)

    # plot
    left = bins[:-1] + start
    bottom = 0
    width = np.diff(bins)
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 1))
    if colors is None:
        colors = 'grey', 'w'
    else:
        assert len(colors) == 2, 'bad colors'
    for h, c, l in zip([h_n, h_a], colors, ['inaccessible', 'accessible']):
        ax.bar(left, h, width=width, bottom=bottom, color=c, align='edge', label=l)
        bottom += h
    ax.set_xlim(start, end)
    ax.set_yticks(ax.get_ylim())
    ax.set_yticklabels(['0%', '100%'])
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in ax.get_xticks()])

    # legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:-1], labels[:-1], loc='center left', bbox_to_anchor=(1, .5), borderpad=0)



# %% [markdown]
# ## Plot repeats

# %%
repeats = geneset_to_pandas(allel.FeatureTable.from_gff3(os.path.join(phase2_ar1.genome_dir, 'agamP4',
                                                                      'Anopheles-gambiae-PEST_REPEATFEATURES_AgamP4.sorted.gff3.gz'),
                                                         attributes=['Name', 'class', 'repeat_consensus', 'Alias', 'ID']))
repeats.head()

# %%
seqid, start, end = region_vgsc
start = start - 5000
end = end + 5000
repeats_vgsc = repeats.query("(seqid == %r) & (end >= %s) & (start <= %s)" % (seqid, start, end)).copy()
repeats_vgsc.head()

# %%
repeats_vgsc['class'].value_counts()


# %%
def plot_repeats(features, y=0, height=1, label=None, ax=None, color=None, color_key=None, skip_label=None, label_rotation=90, yt=.8):
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 1))
        
    n_features = len(features)
    
    for i, (_, feature) in enumerate(features.iterrows()):
        
        # feature coords
        x = feature.start
        width = feature.end - feature.start + 1
        if color:
            if callable(color_key):
                k = color_key(feature)
            else:
                k = feature[color_key]
            try:
                c = color[k]
            except KeyError:
                c = 'grey'
        else:
            c = 'grey'
        patch = plt.Rectangle((x, y), width, height, lw=.2, edgecolor=c, facecolor=c)
        ax.add_patch(patch)

        if label:
            # text coords
            xt = i / n_features
            s = feature[label]
            if skip_label and s in skip_label:
                pass
            else:
                ax.annotate(s, xy=(x + width/2, y + height), xycoords='data', xytext=(xt, yt), textcoords='axes fraction', ha='center', va='bottom', 
                            arrowprops=dict(arrowstyle='-', shrinkA=1, shrinkB=0, lw=.5, color='#aaaaaa', relpos=(0.5, 0)), rotation=label_rotation)

    ax.set_yticks([])
    return ax



# %%
def _plot():

    padding = 5000
    fig, ax = plt.subplots(figsize=(15, 1))
    sns.despine(ax=ax, left=True)
    palette = sns.color_palette('Set2', 9)
    repeat_color = {
        'LINE': palette[0],
        'SINE?': palette[1],
        'LTR': palette[2],
        'DNA': palette[3],
        'RC': palette[3],
        'MITE': palette[3],

    }
    plot_repeats(repeats_vgsc[repeats_vgsc.source == 'RepeatMasker'], ax=ax, label='class',
                 color=repeat_color, color_key=lambda v: v['class'].split('/')[0])
    ax.set_xlim(region_vgsc.start - padding, region_vgsc.end + padding)
    labels = 'LINE', 'SINE', 'LTR', 'Class II (DNA)'
    handles = [plt.Rectangle((0, 0), 1, 1, color=palette[i]) for i in [0, 1, 2, 3]]
    ax.legend(handles=handles, labels=labels, loc='center left', bbox_to_anchor=(1, .5))
    ax.set_ylim(-.2, 5);
    
_plot()

# %%
repeats_vgsc[repeats_vgsc.source != 'RepeatMasker'].head(10)


# %%
def plot_all_repeats(seqid, start, end, label_tes=True, ax=None, y_rm=2, y_trf=1, y_dust=0, height=.9, skip_label=None, label_rotation=90, yt=.8):

    repeats_sel = repeats.query("(seqid == %r) & (end >= %s) & (start <= %s)" % (seqid, start, end)).copy()

    if ax is None:
        fig, ax = plt.subplots(figsize=(15, 1))
        sns.despine(ax=ax, left=True)

    palette = sns.color_palette('Set2', 9)

    # DUST, TRF
    repeat_color = {
        'dust': palette[4],
        'TRF': palette[5],    
    }
    plot_repeats(repeats_sel[repeats_sel.source == 'dust'], y=y_dust, height=height, ax=ax, 
                 color=repeat_color, color_key='source')
    plot_repeats(repeats_sel[repeats_sel.source == 'TRF'], y=y_trf, height=height, ax=ax, 
                 color=repeat_color, color_key='source')

    # repeat-masker
    repeat_color = {
        'LINE': palette[0],
        'SINE?': palette[1],
        'LTR': palette[2],
        'DNA': palette[3],
        'RC': palette[3],
        'MITE': palette[3],

    }
    if label_tes:
        label = 'class'
    else:
        label = None
    plot_repeats(repeats_sel[repeats_sel.source == 'RepeatMasker'], y=y_rm, height=height, ax=ax, label=label,
                 skip_label=skip_label, label_rotation=label_rotation, yt=yt,
                 color=repeat_color, color_key=lambda v: v['class'].split('/')[0])

    labels = 'LINE', 'SINE', 'LTR', 'Class II (DNA)', 'low complexity (DUST)', 'tandem repeat (TRF)'
    handles = [plt.Rectangle((0, 0), 1, 1, color=palette[i]) for i in [0, 1, 2, 3, 4, 5]]
    ax.legend(handles=handles, labels=labels, loc='center left', bbox_to_anchor=(1, .5), ncol=2)
    ax.set_xlim(start, end)
    if label_tes:
        ax.set_ylim(-.2, 10)
    else:
        ax.set_ylim(-.2, 3)


# %%
plot_all_repeats('2L', region_vgsc.start - 5000, region_vgsc.end + 5000)

# %%
plot_all_repeats('2L', region_vgsc.start - 17000, region_vgsc.end + 17000)

# %%
plot_all_repeats('2L', region_vgsc.start - 5000, region_vgsc.end + 5000, label_tes=False)

# %%
plot_all_repeats('2L', region_vgsc.start - 5000, region_vgsc.end + 5000, label_tes=False, y_rm=0, y_dust=0, y_trf=0, height=1)

# %%
plot_all_repeats('2L', region_vgsc.start - 5000, region_vgsc.start + 5000, label_tes=True, y_rm=0, y_dust=0, y_trf=0, height=1)

# %% [markdown]
# TODO: fix repeat and exon plotting so only annotate if within bounds.

# %% [markdown]
# ## Make a figure

# %%
padding = 5000

genome = phase2_ar1.genome_agamp3
accessibility = phase1_ar3.accessibility
seqid = region_vgsc.seqid
start = region_vgsc.start - padding
end = region_vgsc.end + padding

fig = plt.figure(figsize=(15, 4))
gs = mpl.gridspec.GridSpec(4, 1, height_ratios=[1, .4, .2, .4])
row = -1

# repeats
#########
row += 1
ax = fig.add_subplot(gs[row, 0])
sns.despine(ax=ax, left=True, bottom=True)
plot_all_repeats(seqid, start, end, label_tes=True, ax=ax)
ax.set_xticks([])
ax.set_xlim(start, end)

# exons
#######
row += 1
ax = fig.add_subplot(gs[row, 0])
sns.despine(ax=ax, left=True, bottom=True)
plot_davies_exons(start=start, end=end, ax=ax)
ax.set_xticks([])
ax.set_xlim(start, end)
ax.set_ylim(-.1, 1.4)

# accessibility
###############
row += 1
ax = fig.add_subplot(gs[row, 0])
plot_accessibility(accessibility, seqid, start, end, ax=ax, n_bins=1000)
ax.set_xticks([])
ax.set_xlim(start, end)

# sequence composition
######################
row += 1
ax = fig.add_subplot(gs[row, 0])
plot_seq_composition(genome, seqid, start, end, ax=ax, n_bins=1000)
ax.set_xlim(start, end);


# %%
padding = 10000

genome = phase2_ar1.genome_agamp3
accessibility = phase1_ar3.accessibility
seqid = region_vgsc.seqid
start = region_vgsc.start - padding
end = region_vgsc.end + padding

fig = plt.figure(figsize=(15, 3))
gs = mpl.gridspec.GridSpec(3, 1, height_ratios=[2, .5, 1])
row = -1

# repeats
#########
row += 1
ax = fig.add_subplot(gs[row, 0])
sns.despine(ax=ax, left=True, bottom=True)
plot_davies_exons(start=start, end=end, ax=ax, label_position='bottom')
ax.set_ylim(-.7, 2)
ax = ax.twinx()
sns.despine(ax=ax, left=True, bottom=True)
plot_all_repeats(seqid, start, end, label_tes=True, ax=ax, y_dust=0, y_rm=0, y_trf=0, height=1)
ax.set_xticks([])
ax.set_xlim(start, end)
ax.set_ylim(-.7, 2)

# accessibility
###############
row += 1
ax = fig.add_subplot(gs[row, 0])
plot_accessibility(accessibility, seqid, start, end, ax=ax, n_bins=1000)
ax.set_xticks([])
ax.set_xlim(start, end)

# sequence composition
######################
row += 1
ax = fig.add_subplot(gs[row, 0])
plot_seq_composition(genome, seqid, start, end, ax=ax, n_bins=1000)
ax.set_xlim(start, end);


# %% [markdown]
# ## Setup alignment stats

# %%
def build_zarr_alignment_stats(sample_id, chrom, stats_type, field, shape, dtype, **kwargs):

    grp = zarr.open_group('/kwiat/vector/ag1000g/release/phase2.AR1/alignment_stats/stats.zarr2', mode='a')
    arr = grp.require_dataset('%s/%s' % (chrom, field), shape=shape, dtype=dtype, **kwargs)
    
    # check to see if already built...
    try:
        check = arr.attrs[sample_id]
        log(sample_id, 'skipping', chrom, stats_type, field)
    except KeyError:
        log(sample_id, 'building', chrom, stats_type, field)
        try:
            h5_fn = '/media/aliman/aliman/vector/ag1000g/release/phase2-AR1/alignment_stats/%s.%s.%s.h5' % (sample_id, stats_type, chrom)
            h5f = h5py.File(h5_fn, mode='r')
        except OSError as e:
            log(sample_id, e)
        else:
            x = h5f['data'][field]
            col_idx = sample_ids.index(sample_id)
            arr[:, col_idx] = x
            # mark success
            arr.attrs[sample_id] = True
            log(arr)
    
    

# %%
chrom = '2L'
stats_type = 'mapq'
field = 'reads_all'
shape = len(phase2_ar1.genome_agamp3[chrom]), len(sample_ids)
dtype = np.dtype('i4')
chunk_width = 10
chunk_size = 2**23
chunk_height = int(chunk_size / (chunk_width * dtype.itemsize))
chunks = chunk_height, chunk_width
compressor = zarr.Blosc(cname='zstd', clevel=1, shuffle=False)
order = 'F'

for sample_id in sample_ids:
    build_zarr_alignment_stats(sample_id, chrom, stats_type=stats_type, field=field, shape=shape, dtype=dtype, chunks=chunks, compressor=compressor, order=order)

# %%
alignment_stats = zarr.open_group('/kwiat/vector/ag1000g/release/phase2.AR1/alignment_stats/stats.zarr2', mode='r')
alignment_stats

# %%
a = alignment_stats['2L/reads_all'][:, 0]
a

# %%
cache 


# %%
@cache.memoize
def average_coverage(seqid, sample_id):
    reads_all = alignment_stats[seqid]['reads_all']
    # check to see if built
    if sample_id in reads_all.attrs:
        sample_idx = sample_ids.index(sample_id)
        x = reads_all[:, sample_idx]
        xnz = x[x > 0]
        mean = np.mean(xnz)
        median = np.median(xnz)
        mode = np.argmax(np.bincount(xnz))
        std = np.std(xnz)
        return mean, median, mode, std
    else:
        raise ValueError('stats not built')



# %%
average_coverage('2L', 'AA0040-C')

# %%
selected_sample_ids = []
for s in phase1_ar3.sample_ids:
    try:
        # for now, only use ones that have alignment data
        average_coverage('2L', s)
        selected_sample_ids.append(s)
    except ValueError:
        pass

# %%
# check skew
blacklist_sample_ids = []
v = np.array([average_coverage('2L', s) for s in selected_sample_ids])
fig, ax = plt.subplots(figsize=(3, 3))
ax.plot(v[:, 0], v[:, 2], linestyle=' ', marker='.', mfc='none')
ax.set_xlabel('Mean coverage')
ax.set_ylabel('Modal coverage')
for m, d, s in zip(v[:, 0], v[:, 2], selected_sample_ids):
    if m - d > (.05 * m):
        ax.annotate(s, xy=(m, d), xytext=(10, 0), textcoords='offset points', arrowprops=dict(arrowstyle='-'))
        blacklist_sample_ids.append(s)
ax.set_xlim(left=0)
ax.set_ylim(bottom=0)

# %%
whitelist_sample_ids = [s for s in selected_sample_ids if s not in blacklist_sample_ids]
len(whitelist_sample_ids)

# %%
whitelist_sample_ids[:5]

# %%
whitelist_samples = phase1_ar3.samples.set_index('ox_code').loc[whitelist_sample_ids]
whitelist_samples.head()

# %%
n = 28
ds_sample_ids = []
for pop in phase1_ar3.pop_labels:
    if pop != 'colony':
        l = whitelist_samples[whitelist_samples.population == pop].index.values.tolist()
        ds_sample_ids.extend(random.sample(l, n))
# scramble order
ds_sample_ids = random.sample(ds_sample_ids, len(ds_sample_ids))
len(ds_sample_ids)


# %%
@cache.memoize
def normalized_coverage(seqid, start, end, selected_sample_ids, average='mode', convolve=None):

    # extract raw coverage
    loc = slice(start - 1, end)
    reads_all = alignment_stats[seqid]['reads_all'][loc, :]

    # extract for selected samples
    sample_idxs = [sample_ids.index(s) for s in selected_sample_ids]
    reads_all_sel = np.take(reads_all, sample_idxs, axis=1)

    # normalize 
    norm_idx = ['mean', 'median', 'mode'].index(average)
    norm = np.array([average_coverage(seqid, s)[norm_idx] for s in selected_sample_ids])
    reads_all_norm = reads_all_sel * 2 / norm[None, :]
    
    # smooth
    if convolve is not None:
        out = np.zeros_like(reads_all_norm)
        for i in range(out.shape[1]):
            out[:, i] = np.convolve(reads_all_norm[:, i], convolve, mode='same')
    else:
        out = reads_all_norm

    return out


# %%
x = normalized_coverage(seqid=region_vgsc.seqid, start=region_vgsc.start-5000, end=region_vgsc.end+5000, selected_sample_ids=ds_sample_ids)
x.shape

# %%
x = normalized_coverage(seqid=region_vgsc.seqid, start=region_vgsc.start-5000, end=region_vgsc.end+5000, selected_sample_ids=ds_sample_ids, convolve=np.ones(10)/10)
x.shape

# %%
x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
np.convolve(x, np.ones(3)/3, mode='same')

# %%
samples.population.unique()


# %%
def plot_normalized_coverage(seqid, start, end, selected_sample_ids, average='mode', ax=None, plot_kwargs=None, convolve=None, color_field='population', colors=None):
 
    # setup
    if colors is None:
        colors = phase1_ar3.pop_colors
    samples = phase1_ar3.samples.set_index('ox_code')
    color_keys = samples[color_field][selected_sample_ids]
    reads_all_norm = normalized_coverage(seqid=seqid, start=start, end=end, selected_sample_ids=selected_sample_ids, average=average, convolve=convolve)

    if ax is None:
        fig, ax = plt.subplots(figsize=(15, 4))

    x = np.arange(start, end + 1)
    if plot_kwargs is None:
        plot_kwargs = dict()
    plot_kwargs.setdefault('lw', .1)
    plot_kwargs.setdefault('linestyle', '-')
    for i, (s, k) in enumerate(zip(selected_sample_ids, color_keys)):
        y = reads_all_norm[:, i]
        ax.plot(x, y, color=colors[k], **plot_kwargs)

    ax.set_ylim(0, 12)
    ax.set_yticks(range(13))
    ax.grid(axis='y')
    
    # legend
    labels = color_keys.unique().tolist()
    handles = [plt.Rectangle((0, 0), 1, 1, color=colors[l]) for l in labels] 
    ax.legend(loc='center left', bbox_to_anchor=(1, .5), handles=handles, labels=labels)
    
    return ax


# %%
seqid = region_vgsc.seqid
start = region_vgsc.start - 25000
end = region_vgsc.end + 25000

ax = plot_normalized_coverage(seqid, start, end, ds_sample_ids)
ax.set_xlim(start, end);

# %%
seqid = region_vgsc.seqid
start = region_vgsc.start - 25000
end = region_vgsc.end + 25000

ax = plot_normalized_coverage(seqid, start, end, ds_sample_ids, convolve=np.ones(10)/10)
ax.set_xlim(start, end);

# %%
seqid = region_vgsc.seqid
start = region_vgsc.start - 25000
end = region_vgsc.end + 25000

ax = plot_normalized_coverage(seqid, start, end, ds_sample_ids, convolve=np.ones(100)/100)
ax.set_xlim(start, end);

# %%
seqid = region_vgsc.seqid
start = region_vgsc.start - 25000
end = region_vgsc.end + 25000

ax = plot_normalized_coverage(seqid, start, end, ds_sample_ids, convolve=np.ones(300)/300)
ax.set_xlim(start, end);


# %%
def plot_diagnostics(genome, seqid, start, end, selected_sample_ids, n_bins=1000, label_tes=True, convolve_coverage=None):

    accessibility = phase1_ar3.accessibility

    fig = plt.figure(figsize=(15, 5.5))
    gs = mpl.gridspec.GridSpec(4, 1, height_ratios=[5, 1.5, .5, .7])
    row = -1

    # coverage
    ##########
    row += 1
    ax = fig.add_subplot(gs[row, 0])
    sns.despine(ax=ax, left=True, bottom=True)
    plot_normalized_coverage(seqid, start, end, selected_sample_ids, ax=ax, convolve=convolve_coverage)
    ax.set_xticks([])
    ax.set_xlim(start, end)
    ax.set_ylim(-0.1, 12)

    # repeats and exons
    ###################
    row += 1
    ax = fig.add_subplot(gs[row, 0])
    sns.despine(ax=ax, left=True, bottom=True)
    plot_davies_exons(start=start, end=end, ax=ax, y=0, height=.7, label_position='bottom')
    plot_all_repeats(seqid, start, end, label_tes=label_tes, y_rm=0, y_dust=0, y_trf=0, height=.7, ax=ax, label_rotation=0)
    ax.set_xticks([])
    ax.set_xlim(start, end)
    if label_tes:
        ax.set_ylim(-.3, 1.5)
    else:
        ax.set_ylim(-.5, 1)

    # accessibility
    ###############
    row += 1
    ax = fig.add_subplot(gs[row, 0])
    plot_accessibility(accessibility, seqid, start, end, ax=ax, n_bins=n_bins)
    ax.set_xticks([])
    ax.set_xlim(start, end)

    # sequence composition
    ######################
    row += 1
    ax = fig.add_subplot(gs[row, 0])
    plot_seq_composition(genome, seqid, start, end, ax=ax, n_bins=n_bins)
    ax.set_xlim(start, end)



# %%
plot_diagnostics(genome, region_vgsc.seqid, region_vgsc.start - 25000, region_vgsc.end + 25000, selected_sample_ids=ds_sample_ids, 
                 n_bins=2000, label_tes=False, convolve_coverage=np.ones(500)/500)

# %%
plot_diagnostics(genome, region_vgsc.seqid, region_vgsc.start - 5000, region_vgsc.end + 5000, selected_sample_ids=ds_sample_ids, 
                 n_bins=2000, label_tes=False, convolve_coverage=np.ones(500)/500)

# %%
plot_diagnostics(genome, region_vgsc.seqid, region_vgsc.start - 5000, region_vgsc.start + 5000, selected_sample_ids=ds_sample_ids, 
                 n_bins=200, convolve_coverage=np.ones(50)/50)

# %%
plot_diagnostics(genome, region_vgsc.seqid, region_vgsc.start + 4900, region_vgsc.start + 15500, selected_sample_ids=ds_sample_ids, 
                 n_bins=200, convolve_coverage=np.ones(50)/50)

# %%
plot_diagnostics(genome, region_vgsc.seqid, region_vgsc.start + 15000, region_vgsc.start + 25000, selected_sample_ids=ds_sample_ids, 
                 n_bins=200, convolve_coverage=np.ones(50)/50)

# %%
plot_diagnostics(genome, region_vgsc.seqid, region_vgsc.start + 22800, region_vgsc.start + 32700, selected_sample_ids=ds_sample_ids, 
                 n_bins=200, convolve_coverage=np.ones(50)/50)

# %%
plot_diagnostics(genome, region_vgsc.seqid, region_vgsc.start + 30000, region_vgsc.start + 40000, selected_sample_ids=ds_sample_ids, 
                 n_bins=200, convolve_coverage=np.ones(50)/50)

# %%
plot_diagnostics(genome, region_vgsc.seqid, region_vgsc.start + 40000, region_vgsc.start + 50200, selected_sample_ids=ds_sample_ids, 
                 n_bins=200, convolve_coverage=np.ones(50)/50)

# %%
plot_diagnostics(genome, region_vgsc.seqid, region_vgsc.start + 48800, region_vgsc.start + 59000, selected_sample_ids=ds_sample_ids, 
                 n_bins=200, convolve_coverage=np.ones(50)/50)

# %%
plot_diagnostics(genome, region_vgsc.seqid, region_vgsc.start + 58000, region_vgsc.start + 68000, selected_sample_ids=ds_sample_ids, 
                 n_bins=200, convolve_coverage=np.ones(50)/50)

# %%
plot_diagnostics(genome, region_vgsc.seqid, region_vgsc.start + 68000, region_vgsc.start + 78000, selected_sample_ids=ds_sample_ids, 
                 n_bins=200, convolve_coverage=np.ones(50)/50)

# %% [markdown]
# ## Add MQ0

# %%
h5ls(accessibility['2L'])

# %%
x = accessibility['2L']['high_mq0'][:]

# %%
plt.hist(x, bins=np.arange(1, 770, 2))
plt.autoscale(axis='both', tight=True);


# %%
def plot_accessibility_metric(accessibility, metric, seqid, start, end, ax=None, color=None):
 
    # extract accessibility metric
    y = accessibility[seqid][metric][:]
    ymax = np.max(y)
    y = y[start-1:end] 

    # plot
    x = np.arange(start, end + 1)
    if ax is None:
        fig, ax = plt.subplots(figsize=(15, 1))
    ax.fill_between(x, 0, y, color=color, label=metric)
    ax.set_xlim(start, end)
    ax.set_yticks([0, ymax])
    ax.set_ylim(0, ymax)
    ax.set_yticklabels(['0%', '100%'])
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in ax.get_xticks()])
    ax.legend(loc='center left', bbox_to_anchor=(1, .5))



# %%
plot_accessibility_metric(accessibility, 'high_mq0', '2L', region_vgsc.start-5000, region_vgsc.end+5000, color='orange')

# %%
plot_accessibility_metric(accessibility, 'low_mq', '2L', region_vgsc.start-5000, region_vgsc.end+5000, color='orange')

# %%
plot_accessibility_metric(accessibility, 'low_coverage', '2L', region_vgsc.start-5000, region_vgsc.end+5000, color='blue')

# %%
plot_accessibility_metric(accessibility, 'high_coverage', '2L', region_vgsc.start-5000, region_vgsc.end+5000, color='red')


# %%
def plot_diagnostics_v2(genome, seqid, start, end, selected_sample_ids, n_bins=1000, label_tes=True, convolve_coverage=None, coverage_color_field='population', coverage_colors=None):

    accessibility = phase1_ar3.accessibility

    fig = plt.figure(figsize=(15, 7))
    gs = mpl.gridspec.GridSpec(8, 1, height_ratios=[5, 1.5, .5, .5, .5, .5, .5, .7], hspace=0.1)
    row = -1

    # coverage
    ##########
    row += 1
    ax = fig.add_subplot(gs[row, 0])
    sns.despine(ax=ax, left=True, bottom=True)
    plot_normalized_coverage(seqid, start, end, selected_sample_ids, ax=ax, convolve=convolve_coverage, color_field=coverage_color_field, colors=coverage_colors)
    ax.set_xticks([])
    ax.set_xlim(start, end)
    ax.set_ylim(-0.1, 12)
    ax.set_ylabel('Normalized coverage')

    # repeats and exons
    ###################
    row += 1
    ax = fig.add_subplot(gs[row, 0])
    sns.despine(ax=ax, left=True, bottom=True)
    plot_davies_exons(start=start, end=end, ax=ax, y=0, height=.7, label_position='bottom')
    plot_all_repeats(seqid, start, end, label_tes=label_tes, y_rm=0, y_dust=0, y_trf=0, height=.7, ax=ax, label_rotation=0, yt=.9)
    ax.set_xticks([])
    ax.set_xlim(start, end)
    if label_tes:
        ax.set_ylim(-.8, 1.5)
    else:
        ax.set_ylim(-.8, 1)

    # accessibility
    ###############

    row += 1
    ax = fig.add_subplot(gs[row, 0])
    plot_accessibility(accessibility, seqid, start, end, ax=ax, n_bins=n_bins)
    ax.set_xticks([])
    ax.set_xlim(start, end)

    row += 1
    ax = fig.add_subplot(gs[row, 0])
    sns.despine(ax=ax)
    plot_accessibility_metric(accessibility, 'high_mq0', seqid, start, end, ax=ax, color='orange')
    ax.set_xticks([])
    ax.set_xlim(start, end)
    ax.axhline(y=1, linestyle=':', color='grey')
    ax.set_ylim(0, 10)
    ax.set_yticks([0, 10])
    ax.set_yticklabels([0, 10])

    row += 1
    ax = fig.add_subplot(gs[row, 0])
    sns.despine(ax=ax)
    plot_accessibility_metric(accessibility, 'no_coverage', seqid, start, end, ax=ax, color='#333333')
    ax.set_xticks([])
    ax.set_xlim(start, end)
    ax.axhline(y=1, linestyle=':', color='grey')
    ax.set_ylim(0, 10)
    ax.set_yticks([0, 10])
    ax.set_yticklabels([0, 10])

    row += 1
    ax = fig.add_subplot(gs[row, 0])
    sns.despine(ax=ax)
    plot_accessibility_metric(accessibility, 'low_coverage', seqid, start, end, ax=ax, color='blue')
    ax.set_xticks([])
    ax.set_xlim(start, end)
    ax.axhline(y=76, linestyle=':', color='grey')
    ax.set_ylim(0, 200)
    ax.set_yticks([0, 200])
    ax.set_yticklabels([0, 200])

    row += 1
    ax = fig.add_subplot(gs[row, 0])
    sns.despine(ax=ax)
    plot_accessibility_metric(accessibility, 'high_coverage', seqid, start, end, ax=ax, color='red')
    ax.set_xticks([])
    ax.set_xlim(start, end)
    ax.axhline(y=15, linestyle=':', color='grey')
    ax.set_ylim(0, 50)
    ax.set_yticks([0, 50])
    ax.set_yticklabels([0, 50])

    # sequence composition
    ######################
    row += 1
    ax = fig.add_subplot(gs[row, 0])
    plot_seq_composition(genome, seqid, start, end, ax=ax, n_bins=n_bins)
    ax.set_xlim(start, end)
    
    plt.show()



# %%
plot_diagnostics_v2(genome, region_vgsc.seqid, region_vgsc.start - 25000, region_vgsc.end + 25000, selected_sample_ids=ds_sample_ids, 
                 n_bins=2000, label_tes=False, convolve_coverage=np.ones(500)/500)

# %%
plot_diagnostics_v2(genome, region_vgsc.seqid, region_vgsc.start - 5000, region_vgsc.end + 5000, selected_sample_ids=ds_sample_ids, 
                 n_bins=2000, label_tes=False, convolve_coverage=np.ones(500)/500)

# %%
plot_diagnostics_v2(genome, region_vgsc.seqid, region_vgsc.start - 5000, region_vgsc.start + 5000, selected_sample_ids=ds_sample_ids, 
                 n_bins=200, convolve_coverage=np.ones(50)/50)

# %%
plot_diagnostics_v2(genome, region_vgsc.seqid, region_vgsc.start + 4900, region_vgsc.start + 15500, selected_sample_ids=ds_sample_ids, 
                 n_bins=200, convolve_coverage=np.ones(50)/50)

# %%
plot_diagnostics_v2(genome, region_vgsc.seqid, region_vgsc.start + 15000, region_vgsc.start + 25000, selected_sample_ids=ds_sample_ids, 
                 n_bins=200, convolve_coverage=np.ones(50)/50)

# %%
plot_diagnostics_v2(genome, region_vgsc.seqid, region_vgsc.start + 22800, region_vgsc.start + 32700, selected_sample_ids=ds_sample_ids, 
                 n_bins=200, convolve_coverage=np.ones(50)/50)

# %%
plot_diagnostics_v2(genome, region_vgsc.seqid, region_vgsc.start + 30000, region_vgsc.start + 40000, selected_sample_ids=ds_sample_ids, 
                 n_bins=200, convolve_coverage=np.ones(50)/50)

# %%
plot_diagnostics_v2(genome, region_vgsc.seqid, region_vgsc.start + 40000, region_vgsc.start + 50200, selected_sample_ids=ds_sample_ids, 
                 n_bins=200, convolve_coverage=np.ones(50)/50)

# %%
plot_diagnostics_v2(genome, region_vgsc.seqid, region_vgsc.start + 48800, region_vgsc.start + 59000, selected_sample_ids=ds_sample_ids, 
                 n_bins=200, convolve_coverage=np.ones(50)/50)

# %%
plot_diagnostics_v2(genome, region_vgsc.seqid, region_vgsc.start + 58000, region_vgsc.start + 68000, selected_sample_ids=ds_sample_ids, 
                 n_bins=200, convolve_coverage=np.ones(50)/50)

# %%
plot_diagnostics_v2(genome, region_vgsc.seqid, region_vgsc.start + 68000, region_vgsc.start + 78000, selected_sample_ids=ds_sample_ids, 
                 n_bins=200, convolve_coverage=np.ones(50)/50)

# %% [markdown]
# ## Color by kdr

# %%
phase1_ar3.samples.columns

# %%
samples = phase1_ar3.samples[['ox_code', 'population', 'kdr_1014']].set_index('ox_code')
samples.head()

# %%
samples_kdr_ff = samples[samples.kdr_1014 == 'F/F'].index.values.tolist()
len(samples_kdr_ff)

# %%
samples_kdr_ss = samples[samples.kdr_1014 == 'S/S'].index.values.tolist()
len(samples_kdr_ss)

# %%
samples_kdr_wt = samples[samples.kdr_1014 == '+/+'].index.values.tolist()
len(samples_kdr_wt)

# %%
len(whitelist_sample_ids)

# %%
kdr_sample_ids = []
kdr_sample_ids.extend([s for s in samples_kdr_ff if s in whitelist_sample_ids])
kdr_sample_ids.extend([s for s in samples_kdr_ss if s in whitelist_sample_ids])
kdr_sample_ids.extend([s for s in samples_kdr_wt if s in whitelist_sample_ids])
# scramble order
kdr_sample_ids = random.sample(kdr_sample_ids, len(kdr_sample_ids))
len(kdr_sample_ids)

# %%
seqid = region_vgsc.seqid
start = region_vgsc.start - 5000
end = region_vgsc.end + 5000

ax = plot_normalized_coverage(seqid, start, end, kdr_sample_ids, convolve=np.ones(300)/300, color_field='kdr_1014', 
                              colors={'F/F': 'red', 'S/S': 'blue', '+/+': 'green'})
ax.set_xlim(start, end);

# %%
kdr_color_field='kdr_1014'
kdr_colors={'F/F': 'red', 'S/S': 'blue', '+/+': 'green'}

plot_diagnostics_v2(genome, region_vgsc.seqid, region_vgsc.start - 25000, region_vgsc.end + 25000, selected_sample_ids=kdr_sample_ids, 
                 n_bins=2000, label_tes=False, convolve_coverage=np.ones(500)/500, coverage_color_field=kdr_color_field, coverage_colors=kdr_colors)

# %%
plot_diagnostics_v2(genome, region_vgsc.seqid, region_vgsc.start - 5000, region_vgsc.end + 5000, selected_sample_ids=kdr_sample_ids, 
                 n_bins=2000, label_tes=False, convolve_coverage=np.ones(500)/500, coverage_color_field=kdr_color_field, coverage_colors=kdr_colors)

# %%
plot_diagnostics_v2(genome, region_vgsc.seqid, region_vgsc.start - 5000, region_vgsc.start + 5000, selected_sample_ids=kdr_sample_ids, 
                 n_bins=200, convolve_coverage=np.ones(50)/50, coverage_color_field=kdr_color_field, coverage_colors=kdr_colors)

# %%
plot_diagnostics_v2(genome, region_vgsc.seqid, region_vgsc.start + 4900, region_vgsc.start + 15500, selected_sample_ids=kdr_sample_ids, 
                 n_bins=200, convolve_coverage=np.ones(50)/50, coverage_color_field=kdr_color_field, coverage_colors=kdr_colors)

# %%
plot_diagnostics_v2(genome, region_vgsc.seqid, region_vgsc.start + 15000, region_vgsc.start + 25000, selected_sample_ids=kdr_sample_ids, 
                 n_bins=200, convolve_coverage=np.ones(50)/50, coverage_color_field=kdr_color_field, coverage_colors=kdr_colors)

# %%
plot_diagnostics_v2(genome, region_vgsc.seqid, region_vgsc.start + 22800, region_vgsc.start + 32700, selected_sample_ids=kdr_sample_ids, 
                 n_bins=200, convolve_coverage=np.ones(50)/50, coverage_color_field=kdr_color_field, coverage_colors=kdr_colors)

# %%
plot_diagnostics_v2(genome, region_vgsc.seqid, region_vgsc.start + 30000, region_vgsc.start + 40000, selected_sample_ids=kdr_sample_ids, 
                 n_bins=200, convolve_coverage=np.ones(50)/50, coverage_color_field=kdr_color_field, coverage_colors=kdr_colors)

# %%
plot_diagnostics_v2(genome, region_vgsc.seqid, region_vgsc.start + 35000, region_vgsc.start + 45300, selected_sample_ids=kdr_sample_ids, 
                 n_bins=200, convolve_coverage=np.ones(50)/50, coverage_color_field=kdr_color_field, coverage_colors=kdr_colors)

# %%
plot_diagnostics_v2(genome, region_vgsc.seqid, region_vgsc.start + 40000, region_vgsc.start + 50200, selected_sample_ids=kdr_sample_ids, 
                 n_bins=200, convolve_coverage=np.ones(50)/50, coverage_color_field=kdr_color_field, coverage_colors=kdr_colors)

# %%
plot_diagnostics_v2(genome, region_vgsc.seqid, region_vgsc.start + 48800, region_vgsc.start + 59000, selected_sample_ids=kdr_sample_ids, 
                 n_bins=200, convolve_coverage=np.ones(50)/50, coverage_color_field=kdr_color_field, coverage_colors=kdr_colors)

# %%
plot_diagnostics_v2(genome, region_vgsc.seqid, region_vgsc.start + 58000, region_vgsc.start + 68000, selected_sample_ids=kdr_sample_ids, 
                 n_bins=200, convolve_coverage=np.ones(50)/50, coverage_color_field=kdr_color_field, coverage_colors=kdr_colors)

# %%
plot_diagnostics_v2(genome, region_vgsc.seqid, region_vgsc.start + 68000, region_vgsc.start + 78000, selected_sample_ids=kdr_sample_ids, 
                 n_bins=200, convolve_coverage=np.ones(50)/50, coverage_color_field=kdr_color_field, coverage_colors=kdr_colors)

# %%
