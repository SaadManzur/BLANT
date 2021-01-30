import argparse

from utils.read_file import read_file
from utils.score import get_auroc
from utils.score import get_aupr
from utils.score import get_ndcg

def get_parser():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--has-header', action='store_true')
    parser.add_argument('--fpath', type=str, help="Path to the tsv file")
    parser.add_argument('--out', type=str, help="Output path prefix for the images")
    
    return parser

if __name__ == "__main__":
    
    parser = get_parser()
    
    args = parser.parse_args()
    
    if not args.fpath and not args.out:
        
        print("* You either did not specify input file path or output prefix")
        parser.print_help()
    
    else:
        result = read_file(args.fpath, args.has_header)

        print(f"AUROC: {get_auroc(result, args.out)}")

        print(f"AUPR: {get_aupr(result, args.out)}")

        print(f"NDCG: {get_ndcg(result)}")
    