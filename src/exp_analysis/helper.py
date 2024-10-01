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
def calculate_feature_distance(sample: pd.DataFrame, control: pd.DataFrame, top_n=None):
    ''' calculates eclidean distance for given features'''
    if top_n is not None:
        sample = sample.sort_values(by='feature_rank')[:top_n] # to control the noise 
    
    entries_common = pd.Index(np.intersect1d(sample.index, control.index)) 
    control = control.reindex(entries_common)
    # - for feature (raw), missing values are 0
    control['feature'] = control['feature'].fillna(0)
    # - for feaure rank, missing values are maximum
    control['feature_rank'] = control['feature_rank'].fillna(max(control['feature_rank']))
    
    distance_raw = sample['feature'] - control['feature'] # euclidean  
    distance_rank = control['feature_rank'] -sample['feature_rank'] # euclidean  
    
    distance_df = pd.DataFrame({'distance_raw':distance_raw, 'distance_rank':distance_rank}, index=entries_common)
    return distance_df 
def cosine_similarity(nets_dict, col_name='source', weight_col='weight', figsize=(4, 4)):
    from itertools import combinations
    from sklearn.metrics.pairwise import cosine_similarity
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np
    
    # 1. Extract the source-target-weight triples from each network
    nets_names = list(nets_dict.keys())
    nets = list(nets_dict.values())

    # 2. Create a union of all possible links (source-target) to form a "vector space"
    all_links = set()
    for net in nets:
        all_links.update(net[col_name])
    all_links = sorted(list(all_links))  # Sort to keep consistent ordering

    # 3. Create a matrix where each entry represents the weight of a link for each network
    weighted_matrix = np.zeros((len(nets), len(all_links)))
    
    for i, net in enumerate(nets):
        net_links = net[col_name]
        link_to_weight = {link: weight for link, weight in zip(net_links, net[weight_col])}
        for j, link in enumerate(all_links):
            weighted_matrix[i, j] = link_to_weight.get(link, 0)  # Use 0 if link not present in net

    # 4. Compute the pairwise cosine similarity between all networks based on weights
    cosine_sim_matrix = cosine_similarity(weighted_matrix)
    for i in range(cosine_sim_matrix.shape[0]):
        for j in range(cosine_sim_matrix.shape[1]):
            if i>=j:
                cosine_sim_matrix[i,j]=np.NaN

    # 5. Visualize the Cosine Similarity matrix as a heatmap
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    
    sns.heatmap(cosine_sim_matrix, annot=True, cmap="coolwarm", xticklabels=nets_names, yticklabels=nets_names, ax=ax)
    ax.grid(True)
    ax.set_title('Cosine Similarity')

    # Rotate x labels for readability
    plt.xticks(rotation=45, ha='right')
    
    return cosine_sim_matrix, fig

def jaccard_similarity(nets_dict, col_name='link', figsize=(4, 4)):
    from itertools import combinations
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np
    # 1. Extract the source-target pairs from each network as a set of tuples
    nets_names = list(nets_dict.keys())
    nets = list(nets_dict.values())
    if col_name is None:
        link_sets = [set(net.index.values) for net in nets]
    else:
        link_sets = [set(net[col_name].unique()) for net in nets]
    # 2. Initialize an empty matrix for storing pairwise Jaccard similarities
    n = len(nets)
    jaccard_matrix = np.zeros((n, n))

    # 3. Compute the pairwise Jaccard similarity
    for i, j in combinations(range(n), 2):
        A = link_sets[i]
        B = link_sets[j]
        intersection = len(A.intersection(B))
        union = len(A.union(B))
        jaccard_similarity = intersection / union if union != 0 else 0
        jaccard_matrix[i, j] = jaccard_similarity
        jaccard_matrix[j, i] = np.NaN
    # Fill diagonal with 1s (as similarity of a network with itself is 1)
    np.fill_diagonal(jaccard_matrix, np.NaN)
    # 4. Visualize the Jaccard matrix as a heatmap
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    sns.heatmap(jaccard_matrix, annot=True, cmap="coolwarm", xticklabels=nets_names, yticklabels=nets_names, ax=ax)
    ax.grid(True)
    ax.set_title('Jaccard Similarity')
    # Rotate x labels for readability
    plt.xticks(rotation=45, ha='right')

    return jaccard_matrix, fig

def plot_cumulative_density(data, title='', ax=None, s=1, x_label='Weight', **kwdgs):
    # Sort the data
    sorted_data = np.sort(data)
    # Compute the cumulative density values
    cdf = np.arange(1, len(sorted_data) + 1) / len(sorted_data)
    if ax is None:
    	fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    else:
    	fig = None
    ax.step(sorted_data, cdf, where='post', label=title, **kwdgs)
    ax.set_xlabel(x_label)
    ax.set_ylabel('Cumulative Density')
    ax.set_title(title)
    ax.grid(True)
    return fig, ax

class Connectivity:
    def __init__(self, net, **kwargs):
        self.out_deg = determine_centrality(net, source='source', target='target',  **kwargs)
        self.in_deg = determine_centrality(net, source='target', target='source',  **kwargs)
        
class VSA_analysis:
    '''Vester's sensitivity analysis. 
                - active_sum: active sum. Sum along rows of the influence matrix and it indicates how much does a variable influence all the others.
                - passive_sum: passive sum. Its is the sum along columns of the influence matrix and it indicates how sensitive a variable is, how does it react to the influence of others
    '''
    def __init__(self, sample: pd.DataFrame, control: pd.DataFrame, mode='weight', top_q_net = 0.75, critical_change_q_t: float = 0.75):
        # - keep only top quantile links
        control = control[control.weight > control.weight.quantile(top_q_net)]
        sample = sample[sample.weight > sample.weight.quantile(top_q_net)]
        
        print('Determine roles')
        self.vsa_control = self.determine_roles(control)
        self.vsa_sample = self.determine_roles(sample)
        print('Determine role change')
        self.oo = self.diff_roles(self.vsa_control, self.vsa_sample, critical_change_q_t)
        
    @staticmethod
    def determine_roles(net, mode='weight', top_q_role=0.9) -> pd.DataFrame:
        # active sum and passive sum 
        gene_names = np.union1d(net.source.unique(), net.target.unique())
        if mode=='weight':
            active_sum = np.asarray([sum(net.query(f"source == '{gene}'")['weight']) for gene in gene_names])
            passive_sum = np.asarray([sum(net.query(f"target == '{gene}'")['weight']) for gene in gene_names])
        else:
            raise ValueError('define first')

        # define the roles ['Buffering', 'Passive', 'Active', 'Critical'] -> [0, 1, 2, 3]
        active_sum_threhsold = np.quantile(active_sum, top_q_role)
        passive_sum_threhsold = np.quantile(passive_sum, top_q_role)
        # active_sum, passive_sum = standardize_array(active_sum), standardize_array(passive_sum)
        roles = [2*as_flag+ps_flag for as_flag, ps_flag in zip(active_sum>active_sum_threhsold, passive_sum>passive_sum_threhsold)]
        roles = pd.DataFrame(data={'gene_name': gene_names, 'active_sum': active_sum, 'passive_sum': passive_sum, 'role':roles })
        roles.set_index('gene_name',inplace=True)
        return roles
    @staticmethod
    def diff_roles(control: pd.DataFrame, sample: pd.DataFrame, critical_change_q_t: float=.75) -> pd.DataFrame:
        '''
        Find the distance in the role from control to sample, and flags the top changes.
        '''
        # - match the index
        common_genes = pd.Index(np.union1d(sample.index, control.index))
        control = control.reindex(common_genes).fillna(0)
        sample = sample.reindex(common_genes).fillna(0)
        # - calculate distance 
        as_distance = (sample['active_sum'] - control['active_sum'])
        ps_distance = (sample['passive_sum'] - control['passive_sum'])
        overall_distance = as_distance**2 + ps_distance**2
        df_distance = pd.DataFrame({'as_distance':as_distance, 'ps_distance':ps_distance, 'overall_distance':overall_distance}, index=common_genes)

        # df_distance = df_distance.assign(passive_sum_c=control['passive_sum'], active_sum_c=control['active_sum'])
        # df_distance = df_distance.assign(passive_sum_s=sample['passive_sum'], active_sum_s=sample['active_sum'])

        # - define x top percentile as the cut-off value to determine critical role change
        # df_distance['critical_change_as'] = df_distance['as_distance'] > df_distance['as_distance'].quantile(critical_change_q_t)
        # df_distance['critical_change_ps'] = df_distance['ps_distance'] > df_distance['ps_distance'].quantile(critical_change_q_t)
        # df_distance['critical_change_overall'] = df_distance['overall_distance'] > df_distance['overall_distance'].quantile(critical_change_q_t)
        return df_distance
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
    @staticmethod
    def plot_centrality_barh(df, title='',ax=None, xlabel='Degree', ylabel='Gene', colors=None):
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
    def top_edges(self, quantile=0.95):
        grn = self.net 
        grn = grn[grn.weight>grn.weight.quantile(quantile)]
        return grn.link.values
    def calculate_basic_stats(self): 
        self.stats = dict({'n_links': self.net.shape[0], 'n_source': self.net.source.nunique(), 
                    'n_target': self.net.target.nunique(), 'ratio_positive_negative':(self.net.weight>=0).sum()/(self.net.shape[0])},
                    )
        return self.stats
    def calculate_centrality(self) -> None:
        '''Calculate network centrality metrics such as in and out degree
        '''
        self.tf_gene = Connectivity(self.net, normalize=False)
        if self.peak_gene_net is None:
            self.peak_gene = None
        else:
            peak_gene_net = self.peak_gene_net
            peak_gene_net.rename(columns={'peak':'source'},inplace=True)
            self.peak_gene = Connectivity(peak_gene_net, normalize=False)
    def plot_centrality(self, values, **kywds):
        return plot_cumulative_density(values, **kywds)
    def plot_grn_cdf(self, ax=None, title=''):
        values_n = self.net.weight/self.net.weight.max()
        return plot_cumulative_density(values_n, ax=ax, title=title)
    @staticmethod
    def plot_cumulative_density(values, **kywds):
        return plot_cumulative_density(values, **kywds)
    @staticmethod
    def subset_quantile(df, col_name='weight', top_q=0.95, top_n=None, ascending=False):
        if top_n is None:
            df = df[df[col_name] > df[col_name].quantile(top_q)]
        else:
            df = df.sort_values(by=col_name, ascending=ascending, key=abs)[:top_n]
        return df
    def annotate_peaks(self, annotation_df) -> dict[str, float]:
        '''Annotate peaks with associated regions on genome.
        '''
        if self.peak_gene_net is None:
            print('Peak gene net is not given. Peak annoation is skipped.')
            return
        peaks = self.peak_gene_net.peak.unique()
        peaks = self.format_peak(peaks)
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
    @staticmethod
    def calculate_feature(net, feature='active_sum_weight'):
        '''calculate the feature value, e.v. active sum, using net.'''
        if feature=='active_sum_weight':
            net.loc[:, 'weight'] = net.weight.abs()
            df = net.groupby('source')['weight'].sum().to_frame()
        elif feature=='active_sum_degree':
            df = net.groupby('source').size().to_frame()
        elif feature=='passive_sum_degree':
            df = net.groupby('target').size().to_frame()
        else:
            raise ValueError('define first')
        df.columns = ['feature']
        df['feature_rank'] = df.rank(ascending=False)
        return df
    

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
