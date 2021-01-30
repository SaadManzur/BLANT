import argparse

from utils.read_file import read_file
from utils.score import get_auroc
from utils.score import get_aupr
from utils.score import get_ndcg

def get_args():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--has-header', action='store_true')
    parser.add_argument('--fpath', type=str, help="Path to the tsv file")
    parser.add_argument('--out', type=str, help="Output path prefix for the images")
    
    return parser.parse_args()

if __name__ == "__main__":
    
    args = get_args()
    
    result = read_file(args.fpath, args.has_header)
    
    print(f"AUROC: {get_auroc(result, args.out)}")
    
    print(f"AUPR: {get_aupr(result, args.out)}")
    
    print(f"NDCG: {get_ndcg(result)}")
    