import pandas as pd

class infotable(object):
    def __init__(self, target_taxa): 
        self.colnames = ["contig","prot_len","coverage","gc","pid","sp_os","lineage","evalue","parse_lineage"]
        empty_infotable = {}
        for col in self.colnames:
            empty_infotable[col] = []
        
        self.df = pd.DataFrame(empty_infotable, columns=self.colnames)
        if target_taxa is not None:
            self.tar = target_taxa['target']
            self.ex = target_taxa['exception']
        self.keep = None
        self.dump = None
    
    def populate(self, ldict, colnames = None):
        self.df = pd.DataFrame(ldict, columns=colnames)

    def load(self, tsv):
        self.df = pd.read_csv(tsv, sep="\t", header=None)
        self.df.columns = self.colnames

    def clear_decisions(self):
        self.keep = None
        self.dump = None

    def retarget(self, new_target_taxa):
        self.tar = new_target_taxa['target']
        self.ex = new_target_taxa['exception']
    
    def parse_lineage(self):
        self.df = self.df.apply(it_parse_lin, axis=1, args=(self.tar, self.ex))

    def decide_taxonomy (self):
        self.keep = []
        self.dump = []
        grouped = self.df.groupby('contig')
        reformed = grouped.agg({'pid': lambda x: ','.join(x),
                            'evalue': lambda x: ','.join(str(e) for e in x),
                            'parse_lineage': lambda x: ','.join(x)})
    
        reformed.parse_lineage = reformed.parse_lineage.apply(str.split,args=',')
        reformed.pid = reformed.pid.apply(str.split,args=',')
        reformed.evalue = reformed.evalue.apply(str.split,args=',')
        reformed = reformed.reset_index()
    
        for i in reformed.itertuples():
            num_t = i.parse_lineage.count('target')
            num_nt = i.parse_lineage.count('nontarget')
            if num_t > num_nt:
                self.keep.append(i.contig)
            elif num_t == num_nt:
                evalues = [float(e) for e in i.evalue]
                maxes = [x for x in evalues if x == min(evalues)] ## best evalue means lowest, ie min()
                if len(maxes) == 1:
                    best_idx = evalues.index(min(evalues))
                    if i.parse_lineage[best_idx] == 'target':
                        self.keep.append(i.contig)
                    else:
                        self.dump.append(i.contig)
                else:
                    self.dump.append(i.contig)
            else:
                self.dump.append(i.contig)

def it_parse_lin(row, tar, ex):
            if row['lineage'] == "Not_in_taxdb":
                row['parse_lineage'] = 'unclassified'
            elif any(i in tar for i in row['lineage']) is True and any(i in ex for i in row['lineage']) is False:
                row['parse_lineage'] = 'target'
            else:
                row['parse_lineage'] = 'nontarget'
            return row

