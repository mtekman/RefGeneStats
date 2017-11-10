#!/usr/bin/env python3

import sys
import pickle
import os.path

inp = sys.argv[1]

chrom_map = {}  # chr → gene → isoform/canonical → ([exon_sizes], [intron_sizes], direct)
pfile = "pickle.map"

def populateMap():

    # check if pickle exists
    if os.path.exists(pfile):
        print("Previous map data exists, loading from pickle...(delete it to refresh)")
        return 0
        
    
    with open(inp, 'r') as fio:

        fio.readline()

        for line in fio:
            chrom, start, stop, geneinfo, direct, rgname, frames = line.split('\t')
            
            if chrom not in chrom_map:
                chrom_map[chrom] = {}

            isogname, ex_in = geneinfo.split('|')

            length = int(stop) - int(start)
            if length<0:
                length = -length

            gname = isogname.split('-ISO')[0]

            if gname not in chrom_map[chrom]:
                chrom_map[chrom][gname] = {}

            if isogname not in chrom_map[chrom][gname]:
                # exon lengths, intron lengths, directions
                chrom_map[chrom][gname][isogname] = ([],[],[])  


            if ex_in[:2] == "Ex":
                # skip splice sites
                if len(ex_in.split('_')) > 1:
                    continue

                chrom_map[chrom][gname][isogname][0].append(length)
            elif ex_in[:2] == "In":
                chrom_map[chrom][gname][isogname][1].append(length)
            else:
                print("Nuts", line)
                exit(-1)

            direct_bool = direct == "-"

            chrom_map[chrom][gname][isogname][2].append(direct_bool)
            print(chrom, end='\r', file=sys.stderr)

    # Dump map to file for fast reload later
    with open(pfile,'wb') as pio:
        pickle.dump(chrom_map, pio)
        pio.close()



        
class StatCounter:
    """
    Collects data on genes per chromosome
    """

    def __init__(self, chrom):
        self.chrom = chrom
        self.exon_stats = {
            "amount" : [],     # amount of exons per gene
            "min_size" : [],   # min,max,avg size of exons in said genes
            "max_size" : [],
            "avg_size" : []
        } 
        self.intr_stats = {
            "amount" : [],    # amount of introns per gene
            "min_size" : [],  # min,max,avg size of introns in said genes
            "max_size" : [],
            "avg_size" : []
        } 
        self.direct_stats = []  # direction of each exon and intron in said genes
       
    def insertGeneExonicData(self, amount, min_size, max_size, avg_size):
        self.exon_stats["amount"].append(amount)
        self.exon_stats["min_size"].append(min_size)
        self.exon_stats["max_size"].append(max_size)
        self.exon_stats["avg_size"].append(avg_size)

    def insertGeneIntronicData(self, amount, min_size, max_size, avg_size):
        self.intr_stats["amount"].append(amount)
        self.intr_stats["min_size"].append(min_size)
        self.intr_stats["max_size"].append(max_size)
        self.intr_stats["avg_size"].append(avg_size)

    def insertDirect(self, direction):
        self.direct_stats.append(direction)
       
    def writeStats(self, fstream):
        print( self.calcStats(), file=fstream)

    @staticmethod
    def average(array):
        return sum(array)/len(array)

    def calcStats(self):       
        # Number of exons: min, max, average
        exon_amounts = (
            min(self.exon_stats["amount"]),
            max(self.exon_stats["amount"]),
            StatCounter.average(self.exon_stats["amount"])
        )
        # Exon sizes: min, max, average min, average max
        exon_sizes = (
            min(self.exon_stats["min_size"]),
            max(self.exon_stats["max_size"]),
            StatCounter.average(self.exon_stats["avg_size"]),
            StatCounter.average(self.exon_stats["min_size"]),
            StatCounter.average(self.exon_stats["max_size"])
        )

        try:
            intr_amounts = (
                min(self.intr_stats["amount"]),
                max(self.intr_stats["amount"]),
                StatCounter.average(self.intr_stats["amount"])
            )
        except ValueError:
            import pdb; pdb.set_trace()
        
        intr_sizes = (
            min(self.intr_stats["min_size"]),                 # minimum size from list of all minimum sizes
            max(self.intr_stats["max_size"]),                 # maximum size from list of all maximum sizes
            StatCounter.average(self.intr_stats["avg_size"]), # average size from list of all averages
            StatCounter.average(self.intr_stats["min_size"]), # average minimum size (from list of minimums)
            StatCounter.average(self.intr_stats["max_size"])  # average maximum size (from list of maximums)
        )
        
        percent_sense = ( 100 * self.direct_stats.count(True) / len(self.direct_stats) )

        # Num exons (min, max, avg), Size exons (min, max, avg, min_avg, max_avg)
        # Num intrs (min, max, avg), Size intrs (min, max, avg, min_avg, max_avg)
        # % sense
        out_str = self.chrom + "\t"
        out_str += ("%d\t%d\t%.0f" + "\t" + "%d\t%d\t%.0f\t%.0f\t%.0f") % (*exon_amounts, *exon_sizes)
        out_str += "\t"
        out_str += ("%d\t%d\t%.0f" + "\t" + "%d\t%d\t%.0f\t%.0f\t%.0f") % (*intr_amounts, *intr_sizes)
        out_str += "\t"
        out_str += "%.2f" % percent_sense

        return out_str
        
        
        

def processMap():
    global chrom_map
    
    if chrom_map == {}:
        with open(pfile, 'rb') as pio:
            chrom_map = pickle.load(pio)
            pio.close()

    spread_with_isos = open("genes.with_isos.tsv", "w")
    spread_canonical = open("genes.canonical.tsv", "w")

    out_header1 = ("chrom", "gene", "transcript", "sense",
                   "#Exons", "ex(min)", "ex(max)", "ex(avg)",
                   "#Intrs", "in(min)", "in(max)", "in(avg)")
    
    print(*out_header1, sep='\t', file=spread_canonical)
    print(*out_header1, sep='\t', file=spread_with_isos)    

    stats_with_isos = open(spread_with_isos.name.split(".tsv")[0] + '.stats.tsv', 'w')
    stats_canonical = open(spread_canonical.name.split(".tsv")[0] + '.stats.tsv', 'w')

    # Num exons (min, max, avg), Size exons (min, max, avg, min_avg, max_avg)
    # Num intrs (min, max, avg), Size intrs (min, max, avg, min_avg, max_avg)
    # % sense
    out_header_stats1 = ("chrom",
                         "-Number","of","Exons-", "-Size","of","all","the","Exons-",
                         "-Number","of","Intrs-","-Size","of","all","the","Intrs-",
                         "%sense")
    out_header_stats2  = ("     ",
                          "-min", "max", "avg-",
                          "-min", "max", "avg", "min_avg", "max_avg-",
                          "-min", "max", "avg-",
                          "-min", "max", "avg", "min_avg", "max_avg-"
                          )
    
    print(*out_header_stats1, sep='\t', file=stats_with_isos)
    print(*out_header_stats1, sep='\t', file=stats_canonical)
    print(*out_header_stats2, sep='\t', file=stats_with_isos)
    print(*out_header_stats2, sep='\t', file=stats_canonical)

    for chrom in chrom_map:
        # stats
        isochrom = StatCounter(chrom)
        canchrom = StatCounter(chrom)
        
        for gene in chrom_map[chrom]:
            for iso in chrom_map[chrom][gene]:               

                exon_lens, intr_lens, directs = chrom_map[chrom][gene][iso]

                if len(directs) != len(exon_lens) + len(intr_lens):
                    print("mismatched lens", chrom, gene, iso, file=sys.stderr)
                    exit(-1)

                ex_minmaxav = (0, None, None, None)
                in_minmaxav = (0, None, None, None)

                # For current gene:
                if exon_lens != []:                   
                    ex_minmaxav = (
                        len(exon_lens),                 # number of exons
                        min(exon_lens), max(exon_lens), # min/max size of exons
                        sum(exon_lens)/len(exon_lens)   # average size of exons
                    )

                if intr_lens != []:
                    in_minmaxav = (
                        len(intr_lens),                 # number of introns
                        min(intr_lens), max(intr_lens), # min/max size of introns
                        sum(intr_lens)/len(intr_lens)   # average size of introns
                    )

                ddd = list(set(directs))                # direction of gene

                # All exons/introns should be in the same direction as all other
                # exons/introns in the same gene.
                if len(ddd) > 1:
                    print("Nuts in direction", chrom,gene,iso, file=sys.stderr)
                    exit(-1)

                direct = "+" if ddd[0] else "-"
                is_iso = len(iso.split("-ISO")) > 1

                if not(is_iso):
                    if ex_minmaxav[0] != 0:
                        canchrom.insertGeneExonicData(*ex_minmaxav)
                        
                    if in_minmaxav[0] != 0:
                        canchrom.insertGeneIntronicData(*in_minmaxav)
                        
                    canchrom.insertDirect(ddd[0])
                    
                    print(chrom, gene, "CANONICAL", direct, *ex_minmaxav, *in_minmaxav,
                          sep='\t', file=spread_canonical)
                    print(chrom, gene, "CANONICAL", direct, *ex_minmaxav, *in_minmaxav,
                          sep='\t', file=spread_with_isos)
                else:
                    if ex_minmaxav[0] != 0:
                        isochrom.insertGeneExonicData(*ex_minmaxav)
                        
                    if in_minmaxav[0] != 0:
                        isochrom.insertGeneIntronicData(*in_minmaxav)
                        
                    isochrom.insertDirect(ddd[0])
                   
                    print(chrom, gene, iso, direct, *ex_minmaxav, *in_minmaxav,
                          sep='\t', file=spread_with_isos)


        # print chrom stats
        isochrom.writeStats(stats_with_isos)
        canchrom.writeStats(stats_canonical)

        
    spread_with_isos.close()
    spread_canonical.close()

populateMap()
processMap()
        

        

        

        
