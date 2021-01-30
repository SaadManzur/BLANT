class CompiledResult:
    
    def __init__(self):
        
        self.Ts = []
        self.Rs = []
        self.Ps = []
        self.node_pairs = []
        self.confidences = []
        self.scores = []
        self.orbit_pairs = []
        
    def add(self, T, R, P, node_pair, confidence, score, orbit_pair):
        
        self.Ts.append(int(T))
        self.Rs.append(int(R))
        self.Ps.append(float(P))
        self.node_pairs.append(node_pair)
        self.confidences.append(float(confidence))
        self.scores.append(float(score))
        self.orbit_pairs.append(orbit_pair)
        
        
    def get_ts(self):
        
        return self.Ts
    
    def get_scores(self):
        
        return self.scores
    
    def get_Rs(self):
        
        return self.Rs
    
    def __len__(self):
        
        return len(self.Ts)
