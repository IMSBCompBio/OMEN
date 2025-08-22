import yaml

with open('config.yml', 'r') as f:
    config = yaml.safe_load(f)

sample_features_path = config['sample_features_path']
neighbor_matrix_path = config['neighbor_matrix_path']
coords_path = config['coords_path']
