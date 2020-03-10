import pandas as pd
import numpy as np
import ast
import sys
from collections import namedtuple

class InfoTable(object):
    def __init__(self, target_taxa = None, ident="HEAD", colnames = None):
        
        if colnames is None:
            self.colnames = ['contig', 'gc', 'coverage', 'pid', 'length', 'sseqid', 'sp_os', 'desc', 'evalue', 'lineage', 'pertinent_taxlvl', 'parse_lineage']
        else:
            self.colnames = colnames

        empty_infotable = {}
        for col in self.colnames:
            empty_infotable[col] = []

        self.df = pd.DataFrame(empty_infotable, columns=self.colnames)
        self.ttd = target_taxa
        self.tar = None
        self.ex = None
        if target_taxa is not None:
            self.tar = target_taxa['target']
            self.ex = target_taxa['exception']
        self.keep = None
        self.dump = None
        self.ident = ident

        self.children = []
        self.parent = None
        self.depth = 0

    def spawn_child(self, ident, nested_frame):
        children = [c for c in self.children]
        if ident in [c.ident for c in children]:
            duplicate_children = [c for c in children if c.ident == ident]
            assert len(duplicate_children) == 1, "Internatl error, runaway children creation"
            return duplicate_children[0]
        else:
            child = InfoTable(self.ttd, ident=ident)
            child.df = nested_frame
            self.children.append(child)
            child.parent = self
            child.depth = child.parent.depth+1
            return child

    def iter_descendants(self, at_head=True):
        if at_head:
            print (self.ident)
        for c in self.children:
            backbone = 0
            if c.depth > 1:
                backbone = 1
            print ("{}{}|----> {}".format(
                    "|"*(backbone),
                    " "*2*(c.depth-1),
                    c.ident))
            if len(c.children) != 0:
                c.iter_descendants(at_head=False)

    def target_filter(self):
        child_frame = self.df.loc[self.df.parse_lineage == "target"]
        child = self.spawn_child("target", child_frame)
        return child

    def collect_unclassifieds(self, nucl):
        ldict = []
        unclassifieds = [s for s in nucl.seqs() if s.shortname not in self.df.contig.values]
        for u in unclassifieds:
            ldict.append(
                {
                    "contig": u.shortname,
                    "gc": u.gc,
                    "coverage": u.coverage,
                    "taxonomy": "Unclassified"
                }
            )
        unclassifieds = InfoTable()
        unclassifieds.populate(ldict, ldict[0].keys())
        return unclassifieds

    def tidy (self, taxlvl_idx):
        self.df.lineage = self.df.lineage.apply(str.replace,args=(';','.'))

        ### Need to fix this crap and deal with this when building taxdb - too late to be doing this nonsense
        self.df.lineage = self.df.lineage.apply(str.replace,args=(", ",'_'))
        self.df.lineage = self.df.lineage.apply(str.replace, args=(',','_'))

        self.df.lineage = self.df.lineage.apply(str.split,args='.')
        self.df.lineage = self.df.lineage.apply(lambda x: [t.strip() for t in x])

        ## Get taxonomy level now so R doesn't freak out later
        self.df['pertinent_taxlvl'] = self.df.apply(it_get_taxonomy_level, args=(taxlvl_idx,), axis=1)['lineage']
        self.df.reset_index()

    def rfilter(self, feature, window, child_label=""):
        child_frame = self.df.loc[(self.df[feature] >= window[0]) & (self.df[feature] <= window[1])]
        child_ident = "{}:{}{}".format(child_label, feature, window)
        child = self.spawn_child(child_ident, child_frame)
        return child

    def rfilter_inplace_from_parent (self, feature, window):
        assert self.parent is not None, "No parent to filter."
        if window[0] > window[1]:
            window = window[::-1]
        self.df = self.parent.df.loc[(self.parent.df[feature] >= window[0]) & (self.parent.df[feature] <= window[1])]

    def reset_from_parent (self):
        assert self.parent is not None, "No parent to reset from."
        self.df = self.parent.df

    def tnt_population (self):
        TntPop = namedtuple("TntPop",["target","nontarget"])
        num_t = float(sum(self.df.parse_lineage == "target"))
        num_nt = float(sum(self.df.parse_lineage == "nontarget"))
        return TntPop(num_t, num_nt)

    def summary_stats (self, feature):
        Sstats = namedtuple("{}".format(feature),["mean","std","min","max","range"])
        f_mean = np.mean(self.df[feature])
        f_std = np.std(self.df[feature])
        f_min = self.df[feature].min()
        f_max = self.df[feature].max()
        f_range = (f_min,f_max)
        return Sstats(f_mean, f_std, f_min, f_max, f_range)

    def populate(self, ldict, colnames = None):
        self.df = pd.DataFrame(ldict, columns=colnames)

    def load(self, tsv):
        self.df = pd.read_csv(tsv, sep="\t", header=None)
        if self.df.shape[1] == len(self.colnames):
            self.df.columns = self.colnames
        elif self.df.shape[1] == len(self.colnames_new):
            self.df.columns = self.colnames_new
        else:
            self.df.columns = self.unparse_colnames

        self.df.lineage = self.df.lineage.apply(ast.literal_eval)

    def clear_decisions(self):
        self.keep = None
        self.dump = None

    def set_target(self, target, exceptions = None):
        self.tar = [t.strip() for t in target.split(',')]
        if exceptions is None:
            self.ex = []
        else:
            self.ex = [e.strip() for e in exceptions.split(',')]

    def taxon_level(self, level):
        return self.df.apply(it_get_taxonomy_level, axis=1, args=(level,))

    def parse_lineage(self):
        self.df = self.df.apply(it_parse_lin, axis=1, args=(self.tar, self.ex,))

    def decide_inclusion (self):
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


def it_get_taxonomy_level(row, level):
    ret = row.reindex( ['contig','lineage','evalue'] )

    if row['lineage'] == "Not_in_taxdb":
        ret['lineage'] = "unclassified"
    elif len(row['lineage']) <= level:
        ret['lineage'] = "unclassified"
    else:
        ret['lineage'] = ret['lineage'][level]

    return ret

def it_parse_lin(row, tar, ex):
    if row['lineage'] == "Not_in_taxdb":
        row['parse_lineage'] = 'unclassified'
    elif any([i in tar for i in row['lineage']]) and not any([i in ex for i in row['lineage']]):
        row['parse_lineage'] = 'target'
    else:
        row['parse_lineage'] = 'nontarget'
    return row

def get_by_idx (row):
    ret = []
    for i in row.maxes:
        ret.append(row.lineage[i])
    return ret

def count_unique (l):
    if len(l) == 1:
        return l[0]
    else:
        counts = {}
        for ele in set(l):
            counts[l.count(ele)] = ele
        best = {c:ele for c,ele in counts.items() if c == max(counts.keys())}
        if len(best.keys()) > 1:
            #logger.critical("Too many best hits...write more code to deal with this.")
            sys.exit(-5)
        else:
            return best[next(iter(best.keys()))]