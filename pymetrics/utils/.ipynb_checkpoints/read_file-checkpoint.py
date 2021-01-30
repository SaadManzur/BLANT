import csv

from utils.compiled_result import CompiledResult

def read_file(fpath, has_header=False):
    
    result = CompiledResult()
    
    with open(fpath, 'r') as infile_:
        
        file_rows = csv.reader(infile_, delimiter="\t")
        
        for r_count, file_row in enumerate(file_rows):
            
            if has_header and r_count == 0:
                continue
                
            T, R, _, P, node_pair = file_row[0].split(" ")
            confidence = float(file_row[1])
            score, orbit_pair = file_row[2].split(" ")
            
            result.add(T, R, P, node_pair, confidence, score, orbit_pair)
            
        infile_.close()
        
    return result