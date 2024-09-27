import pandas as pd
import anndata as ad
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import minmax_scale



colors_blind = [
    '#E69F00',  # Orange
    '#56B4E9',  # Sky Blue
    '#009E73',  # Bluish Green
    '#F0E442',  # Yellow
    '#0072B2',  # Blue
    '#D55E00',  # Vermillion
    '#CC79A7']  # Reddish Purple


def determine_centrality(net, source='source', target='target', normalize=False, mode='degree') -> pd.DataFrame:
    """Calculates degree of centrality for a given net

    Args:
        net (DataFrame): GRN
        source (str, optional)
        target (str, optional)
        normalize (bool, optional): whether to normalize the centrality to the number of target

    Returns:
        centrality: a df with the size of unique source
    """
    if mode=='degree':
        centrality = net.groupby(source)[target].nunique().to_frame().rename(columns={target:'degree'})
        if normalize:
            total_targets = net[target].nunique()
            centrality = centrality/total_targets
    elif mode=='weight':
        net.weight = net.weight.abs()
        centrality = net.groupby(source)['weight'].sum().to_frame()
        
    else:
        raise ValueError('Define first')
    return centrality


def plot_cumulative_density(data, title='', ax=None, s=1, **kwdgs):
    # Sort the data
    sorted_data = np.sort(data)
    # Compute the cumulative density values
    cdf = np.arange(1, len(sorted_data) + 1) / len(sorted_data)
    if ax is None:
    	fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    else:
    	fig = None
    ax.step(sorted_data, cdf, where='post', label=title, **kwdgs)
    ax.set_xlabel('Weight')
    ax.set_ylabel('Cumulative Density')
    ax.set_title(title)
    ax.grid(True)
    return fig, ax

# class Connectivity:
#     def __init__(self, net, **kwargs):
#         self.out_deg = degree_centrality(net, source='source', target='target',  **kwargs)
#         self.in_deg = degree_centrality(net, source='target', target='source',  **kwargs)
        
class Exp_analysis:
    '''
    This class provides functions for explanatory analysis of GRNs
    '''
    def __init__(self, net, peak_gene_net=None):
        self.net = net 
        self.net.weight = minmax_scale(self.net.weight)
        self.net['link'] = self.net['source'].astype(str) + '_' + self.net['target'].astype(str)

        self.peak_gene_net = peak_gene_net
        self.tfs = net.source.unique()
        self.targets = net.target.unique()
        if peak_gene_net is not None:
            self.cres = peak_gene_net.source.unique()
        # check duplicates
        dup_flag = False
        if 'cell_type' in net.columns:
            if net.duplicated(subset=['source','target','cell_type']).any():
                dup_flag = True
        else:
            if net.duplicated(subset=['source','target']).any():
                dup_flag = True
        if dup_flag:
            raise ValueError('The network has duplicated source target combinations.')
        self.peak_annot = None
    def plot_centrality(self, df, title='',ax=None, xlabel='Degree', ylabel='Gene', colors=None):
        if ax==None:
            fig, ax = plt.subplots(figsize=(10, 6))
        df['degree'].plot(kind='barh', color='skyblue', ax=ax)  # Pass ax to the plot method
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(axis='x', linestyle='--', alpha=0.7)
        if colors is not None:
            for label, color in zip(ax.get_yticklabels(), colors):
                label.set_color(color)
    def role_analysis(self, top_q=None) -> pd.DataFrame:
        '''
            Vester's sensitivity analysis. 
                - active_sum: active sum. Sum along rows of the influence matrix and it indicates how much does a variable influence all the others.
                - passive_sum: passive sum. Its is the sum along columns of the influence matrix and it indicates how sensitive a variable is, how does it react to the influence of others
                - Active: +AS, -PS
                - Passive: -AS, -PS
                - Critical: +AS, +PS
                - neutral: -AS, -PS
        ------------------------------------------------------
        inputs: net (DataFrame)
        top_q: top quantile
        output: VSA (DataFrame) -> genename: active_sum, passive_sum, Role
        '''
        # active sum and passive sum
        net = self.net.copy()  
        net = net[net.weight>net.weight.quantile(top_q)]
        gene_names = np.union1d(net.source.unique(), net.target.unique())
        active_sum = np.asarray([sum(net.query(f"source == '{gene}'")['weight']) for gene in gene_names])
        passive_sum = np.asarray([sum(net.query(f"target == '{gene}'")['weight']) for gene in gene_names])

        # define the roles ['Buffering', 'Passive', 'Active', 'Critical'] -> [0, 1, 2, 3]
        active_sum_threhsold = np.quantile(active_sum, 0.9)
        passive_sum_threhsold = np.quantile(passive_sum, 0.9)
        # active_sum, passive_sum = standardize_array(active_sum), standardize_array(passive_sum)
        roles = [2*as_flag+ps_flag for as_flag, ps_flag in zip(active_sum>active_sum_threhsold, passive_sum>passive_sum_threhsold)]
        self.roles = pd.DataFrame(data={'gene_name': gene_names, 'active_sum': active_sum, 'passive_sum': passive_sum, 'role':roles })
        self.roles.set_index('gene_name',inplace=True)
        return self.roles
    
    def top_edges(self, quantile=0.95):
        grn = self.net 
        grn = grn[grn.weight>grn.weight.quantile(quantile)]
        return grn.link.values
    def calculate_basic_stats(self):
        
        self.stats = dict({'n_links': self.net.shape[0], 'n_source': self.net.source.nunique(), 
                    'n_target': self.net.target.nunique(), 'ratio_positive_negative':(self.net.weight>=0).sum()/(self.net.shape[0])},
                    )
        return self.stats
    def determine_source_properties(self):
        c_degree = determine_centrality(self.net, mode='degree')
        c_weight = determine_centrality(self.net, mode='weight')
        
        self.sources = c_weight.merge(c_degree, left_index=True, right_index=True)
    
    # def calculate_centrality_stats(self) -> None:
    #     '''Calculate network centrality metrics such as in and out degree
    #     '''
    #     self.tf_gene = Connectivity(self.net, normalize=False)
    #     if self.peak_gene_net is None:
    #         self.peak_gene = None
    #         self.n_cres = None
    #     else:
    #         self.peak_gene = Connectivity(self.peak_gene_net, normalize=False)
    #         self.n_cres = self.peak_gene_net.source.nunique()
    def plot_grn_cdf(self, ax=None, title=''):
        values_n = self.net.weight/self.net.weight.max()
        return plot_cumulative_density(values_n, ax=ax, title=title)
    def plot_cumulative_density(self, values, **kywds):
        return plot_cumulative_density(values, **kywds)

    def annotate_peaks(self, annotation_df) -> dict[str, float]:
        '''Annotate peaks with associated regions on genome.
        '''
        peaks = self.format_peak(self.cres)
        annotation_df = annotation_df[annotation_df.peak.isin(peaks)]
        value_counts = annotation_df.annotation.value_counts()
        sum_values = value_counts.sum()
        value_ratio = ((value_counts/sum_values)*100).round(1)

        self.peak_annot = value_ratio.to_dict()
        return self.peak_annot
    @staticmethod
    def format_peak(peaks:list[str]) -> list[str]:
        '''Convert any peak data to chr:start-end
        '''
        import re
        formatted_peaks = []
        for peak in peaks:
            chr_, start, end = re.split(r'[:\-_]', peak)
            peak = f"{chr_}:{start}-{end}"

            formatted_peaks.append(peak)
        return formatted_peaks
    def annotate_genes(self, gene_annotation_df) -> dict[str, float]:
        '''Annotates genes'''
        gene_annot = gene_annotation_df[gene_annotation_df.Gene.isin(self.targets)].Transcript_type.value_counts().to_dict()
        # gene_annot.values
        self.gene_annot = {key:round((value*100/len(self.targets)), 1) for key, value in gene_annot.items()}
        return self.gene_annot
def plot_interactions(interaction_df: pd.DataFrame, min_subset_size=None, min_degree=None, color_map=None) -> plt.Figure: 
    """Upset plot of interactions

    Args:
        interaction_df (pd.DataFrame): _description_

    Returns:
        plt.Figure: _description_
    """
    import upsetplot
    import warnings
    warnings.filterwarnings(action='ignore')
    fig = plt.figure()
    out_dict = upsetplot.plot(upsetplot.from_indicators(indicators=lambda a: a==True, data=interaction_df), fig=fig, 
            show_counts=True, 
            show_percentages = '{:.0%}',
            sort_by='cardinality', 
            # min_subset_size =".1%", # min interaction to show
            min_subset_size = min_subset_size, # min interaction to show
            min_degree=min_degree,
            facecolor='grey',
            other_dots_color=.1, 
            shading_color =.01, 
            with_lines = True,
            element_size = 35,
            intersection_plot_elements=5,
            totals_plot_elements = 2
            )
    # Loop through each bar and change its face color
    matrix_ax, shading_ax, totals_ax, intersections_ax = out_dict['matrix'], out_dict['shading'], out_dict['totals'], out_dict['intersections']

    methods_order = [label.get_text() for label in matrix_ax.get_yticklabels()]
    # methods_order.reverse()
    print(methods_order)
    for i_bar, bar in enumerate(totals_ax.patches):
        if color_map is None:
            bar.set_facecolor(colors_blind[i_bar])
        else:
            bar.set_facecolor(color_map[methods_order[i_bar]])
        bar.set_edgecolor('white')

    for bar in intersections_ax.patches:
        bar.set_facecolor('#c49e81')
        bar.set_edgecolor('black')
        bar.set_linewidth(.4)

    for bar, new_color in zip(shading_ax.patches, colors_blind):
        bar.set_facecolor(new_color)
        bar.set_alpha(.1)
        bar.set_edgecolor('black')

    plt.subplots_adjust(wspace=-.4)  # Example adjustment, modify as needed
    return fig
def create_interaction_df(data_dict: list[str]) -> pd.DataFrame:
    """Receives a dict of name:list and returns a df where the index are unique name in the list, the columns are the names of the data_dict, and the values are booleans, showing the occurance.

    Args:
        data_dict (dict[str, list(str)]): dict of list of names

    Returns:
        interaction: _description_
    """
    
    all_cases = list(set(item for items in data_dict.values() for item in items))
    interaction_df = pd.DataFrame(index=all_cases)
    for key, items in data_dict.items():
        interaction_df[key] = interaction_df.index.isin(items)
    return interaction_df

def create_interaction_info(exp_objs_dict: dict[str, Exp_analysis]) -> dict[str, pd.DataFrame]:
    """Create interaction df for various elements of a grns.

    Args:
        exp_objs_dict (dict[str, Exp_analysis]): a dict with names as grn names and values as Exp_analysis object

    Returns:
        dict[str, pd.DataFrame]: _description_
    """
    
    # interaction of links 
    links_dict = {}
    source_dict = {}
    target_dict = {}
    for name, obj in exp_objs_dict.items():
        # links
        grn = obj.net
        grn['link'] = grn['source'] + '_' + grn['target']
        links_dict[name] = grn['link'].unique()
        # source
        source_dict[name] = grn['source'].unique()
        # target
        target_dict[name] = grn['target'].unique()
    interaction_info = {'links':create_interaction_df(links_dict), 
                        'source':create_interaction_df(source_dict),
                        'target':create_interaction_df(target_dict)}
    return interaction_info
