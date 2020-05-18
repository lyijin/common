#!/usr/bin/env python3

"""
> contrast_bismark_covs.py <

Using seaborn, plot scatterplots to contrast methylation levels / coverage
across multiple cov files. The first cov file provided would be treated as
the baseline for all other cov files to be compared to.

Works with gzipped cov files too!
"""
import argparse
import math

import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats
import seaborn as sns

parser = argparse.ArgumentParser(description="""
Using seaborn, plot scatterplots to contrast methylation levels / coverage
across multiple cov files. The first cov file provided would be treated as
the baseline for all other cov files to be compared to.

Works with gzipped cov files too!""")
parser.add_argument('cov_files', metavar='cov_file',
                    type=argparse.FileType('r'), nargs='+',
                    help='cov files produced by `bismark_methylation_extractor`.')
parser.add_argument('-c', '--mincov', metavar='n', type=int, default=10,
                    help='include positions with coverage >= n in ALL FILES.')
parser.add_argument('-l', '--labels', metavar='label', nargs='+', type=str,
                    help='''override default labeling using filenames with
                            provided labels (one per file).''')
parser.add_argument('--coverage', action='store_true',
                    help='plot coverage plots instead of betas.')
parser.add_argument('--all', action='store_true',
                    help='plots each file vs. all other files.')
args = parser.parse_args()

# sanity checks
#   1. check whether there are at least two files provided
assert len(args.cov_files) > 1, \
    f'At least two files have to be provided! (# files now: {len(args.cov_files)})'

#   2. check whether number of labels tally with number of files
if not args.labels:
    args.labels = [x.name for x in args.cov_files]
else:
    assert len(args.labels) == len(args.cov_files), \
        f'# labels != # files: {len(args.labels)} vs. {len(args.cov_files)}'

#   3. check plot type
plot_type = 'coverage' if args.coverage else 'beta'

def pearson_r_text(x, y, **kws):
    # convert representation of p value to scientific if below 0.0001
    # for scientific: round number to nearest negative power of 10
    m, c, r, p_val, se_m = scipy.stats.linregress(x, y)
    if p_val < 0.0001:
        if p_val:
            p_val = 10 ** (math.ceil(math.log10(p_val)))
        else:
            p_val = 10 ** -323
        p_string = f'p < {p_val:.0e}'
    else:
        p_string = f'p = {p_val:.4f}'
    
    corr_string = f'y = {m:.2f}x {"+" if c >= 0 else "-"} {abs(c):.2f}'
    
    ax = plt.gca()
    ax.annotate(f'{corr_string}\nr = {r:.2f}; {p_string}\nSE(m) = {se_m:.3f}',
                xy=(.95, .95), xycoords=ax.transAxes, 
                verticalalignment='top', horizontalalignment='right')

def read_cov_data(cov_filename, suffix=''):
    '''
    When given a Bismark cov file, read it. Uses pandas to automatically handle
    gzipped file.
    '''
    if suffix != '':
        suffix = '_' + suffix
    
    df = pd.read_table(cov_filename, header=0,
                       names=['scaf', 'pos', 'pos2', f'beta{suffix}',
                              f'meth{suffix}', f'unmeth{suffix}'])
    df = df.drop('pos2', axis=1)
    df[f'coverage{suffix}'] = df[f'meth{suffix}'] + df[f'unmeth{suffix}']
    
    return df

# define df used in the pairplot
plot_df = read_cov_data(args.cov_files[0].name, suffix=args.labels[0])

for n in range(1, len(args.cov_files)):
    temp_df = read_cov_data(args.cov_files[n].name, suffix=args.labels[n])
    plot_df = pd.merge(plot_df, temp_df, how='outer', on=['scaf', 'pos'])

# post-process `plot_df` further. fill NAs with 0, keeping in mind that
# this is OK for coverage, but NOT OK for beta. this hack is overall OK
# because beta plots would require a minimum coverage of 1, which means
# that NA betas converted to 0 won't be plotted.
plot_df = plot_df.fillna(0)

if plot_type == 'beta':
    plot_df['min_coverage'] = plot_df.loc[:, plot_df.columns.str.startswith('coverage_')].min(axis=1)
    plot_df = plot_df[plot_df['min_coverage'] >= args.mincov]

# plot the pairplot!
sns.set(style='ticks', font_scale=1)
#fig, ax = plt.subplots()

if args.all:
    g = sns.pairplot(plot_df.loc[:, plot_df.columns.str.startswith(f'{plot_type}_')],
                     kind='reg', diag_kind='kde', height=3,
                     plot_kws={'line_kws': {'color':'#262626', 'alpha': 0.7}, 
                               'scatter_kws': {'s': 5, 'alpha': 0.2,
                               'edgecolor': 'none'}},
                     diag_kws={'shade': True},
                     grid_kws={'layout_pad': 1})
    
    g.map_offdiag(pearson_r_text)

else:
    g = sns.pairplot(plot_df,
                     kind='reg', height=3,
                     x_vars=[f'{plot_type}_{args.labels[0]}'],
                     y_vars=[f'{plot_type}_{x}' for x in args.labels[1:]],
                     plot_kws={'line_kws': {'color':'#262626', 'alpha': 0.7}, 
                               'scatter_kws': {'s': 5, 'alpha': 0.2,
                               'edgecolor': 'none'}},
                     grid_kws={'layout_pad': 1})

    g.map(pearson_r_text)

# a few other plot type-specific settings
if plot_type == 'beta':
    title_text = f'Betas | min. cov {args.mincov} | n = {len(plot_df):,}'
    
    # always plot in the [0, 100] range
    g.set(xlim=[0, 100], ylim=[0, 100])
    g.set(xticks=[0, 20, 40, 60, 80, 100], yticks=[0, 20, 40, 60, 80, 100])
elif plot_type == 'coverage':
    title_text = f'Coverages | n = {len(plot_df):,}'
    
    # maxminval is the lowest max() across the plotted columns--make dots
    # fill up most of the canvas
    maxminval = plot_df.loc[:, plot_df.columns.str.startswith(f'{plot_type}_')].max().min()
    g.set(xlim=(0, maxminval), ylim=(0, maxminval))
    
    # set y ticks to mirror x ticks
    g.set(yticks=g.fig.axes[0].get_xticks())

g.fig.suptitle(title_text, x=0, y=1, ha='left', va='bottom',
               fontsize='medium', fontweight='bold')

sns.despine(offset=5, trim=True)

# save figure
fig = plt.gcf()

# without bbox_inches, the saved figure has truncated axes.
output_filename = 'output.pdf'
fig.savefig(output_filename, bbox_inches='tight')
