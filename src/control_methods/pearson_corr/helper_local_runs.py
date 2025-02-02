
import argparse
parser = argparse.ArgumentParser(description="Process multiomics RNA data.")
parser.add_argument('--rna', type=str, help='Path to the multiomics RNA file')
parser.add_argument('--prediction', type=str, help='Path to the prediction file')
parser.add_argument('--tf_all', type=str, help='Path to the tf_all')
parser.add_argument('--num_workers', type=str, help='Number of cores')
parser.add_argument('--max_n_links', type=str, help='Number of top links to retain')

args = parser.parse_args()

if args.max_n_links:
    par['max_n_links'] = int(args.max_n_links)
    
if args.rna:
    par['rna'] = args.rna
if args.prediction:
    par['prediction'] = args.prediction
if args.tf_all:
    par['tf_all'] = args.tf_all
if args.num_workers:
    par['num_workers'] = args.num_workers

os.makedirs(par['temp_dir'], exist_ok=True)

meta = {
    'resources_dir': 'src/utils'
    }
sys.path.append(meta["resources_dir"])