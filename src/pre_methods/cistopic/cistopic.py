import pickle
cistopic_obj = pickle.load(
            open(f'../perturb-multiomics-grn/output/infer/scenicplus/scATAC/cistopic_obj.pkl', 'rb'))
# get cell topic association 
cell_topic = cistopic_obj.selected_model.cell_topic.T
cell_names = cistopic_obj.cell_data.obs_id.values

cell_topic.index = cell_names


cell_topic.to_csv(f'resources/grn-benchmark/supp/cell_topic.csv')