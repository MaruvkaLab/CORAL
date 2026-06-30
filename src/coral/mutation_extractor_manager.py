import os
import gzip
import json
import csv
from collections import defaultdict

from .utils import log

REMOVE_CHARS = str.maketrans('', '', '^$[]')
CHR_IDX, POSITION_IDX, REF_NUC_IDX, N_READS_1_IDX, NUC_1_IDX, N_READS_2_IDX, NUC_2_IDX = range(7)
PREV_IDX, CUR_IDX, NEXT_IDX = 0, 1, 2
REF_IDX, TAXA1_IDX, TAXA2_IDX = 0, 1, 2


class MutationExtractor:
    def __init__(self, reference, taxon1, taxon2, pileup_file, mutation_output_dir, triplet_output_dir,
                 no_full_mutations=False, no_cache=False, verbose=True):
        self.reference = reference
        self.taxon1 = taxon1
        self.taxon2 = taxon2
        self.pileup_file = pileup_file
        self.mutation_output_dir = mutation_output_dir
        self.triplet_output_dir = triplet_output_dir
        self.no_full_mutations = no_full_mutations
        self.no_cache = no_cache
        self.verbose = verbose

        self.pileup_file = pileup_file
        self.out_json1 = os.path.join(self.mutation_output_dir, f"{taxon1}__{taxon2}__{reference}__mutations.json")
        self.out_json2 = os.path.join(self.mutation_output_dir, f"{taxon2}__{taxon1}__{reference}__mutations.json")
        self.trip_out_json1 = os.path.join(self.triplet_output_dir, f"{self.taxon1}__{self.taxon2}__{self.reference}__triplets.json")
        self.trip_out_json2 = os.path.join(self.triplet_output_dir, f"{self.taxon2}__{self.taxon1}__{self.reference}__triplets.json")

        self.csv_path1 = None if no_full_mutations else os.path.join(self.mutation_output_dir, f"{taxon1}__{taxon2}__{reference}__mutations.csv.gz")
        self.csv_path2 = None if no_full_mutations else os.path.join(self.mutation_output_dir, f"{taxon2}__{taxon1}__{reference}__mutations.csv.gz")

    def extract(self):
        os.makedirs(self.mutation_output_dir, exist_ok=True)
        os.makedirs(self.triplet_output_dir, exist_ok=True)

        jsons_exist = all(os.path.exists(p) for p in [self.out_json1, self.out_json2, self.trip_out_json1, self.trip_out_json2])
        csvs_exist = (self.no_full_mutations or
                    all(os.path.exists(p) for p in [self.csv_path1, self.csv_path2]))

        if not self.no_cache and jsons_exist and csvs_exist:
            log("Mutation counts already exist. Skipping.", self.verbose)
            return

        species_mut1 = defaultdict(int)
        species_mut2 = defaultdict(int)
        species_triplet1 = defaultdict(int)
        species_triplet2 = defaultdict(int)

        csv1 = csv2 = None
        if not self.no_full_mutations:
            header = "chromosome,position,mutation\n"
            csv1 = gzip.open(self.csv_path1, 'wt')
            csv1.write(header)

            csv2 = gzip.open(self.csv_path2, 'wt')
            csv2.write(header)

        with gzip.open(self.pileup_file, 'rt') as f:

            line_fields = [None, self.parse_line(f.readline()), self.parse_line(f.readline())]
            qc_flags = [False, self.quality_check(line_fields[1]), self.quality_check(line_fields[2])]

            for line in f:
                line_fields = [line_fields[1], line_fields[2], self.parse_line(line)]
                qc_flags = [qc_flags[1], qc_flags[2], self.quality_check(line_fields[2])]

                if all(qc_flags):
                    triplets = self.extract_triplets(line_fields)
                    mut1, mut2, trip1, trip2 = self.detect_mutation_triplet(triplets)
                    chrom = line_fields[1][CHR_IDX]
                    pos = int(line_fields[1][POSITION_IDX])

                    if mut1:
                        species_mut1[mut1] += 1
                        csv1.write(f"{chrom},{pos},{mut1}\n")

                    if mut2:
                        species_mut2[mut2] += 1
                        csv2.write(f"{chrom},{pos},{mut2}\n")

                    if trip1:
                        species_triplet1[trip1] += 1
                    if trip2:
                        species_triplet2[trip2] += 1
        
        if csv1:
            csv1.close()
        if csv2:
            csv2.close()

        with open(self.out_json1, 'w') as f:
            json.dump(species_mut1, f, indent=2)
        with open(self.out_json2, 'w') as f:
            json.dump(species_mut2, f, indent=2)
        
        with open(self.trip_out_json1, 'w') as f:
            json.dump(species_triplet1, f, indent=2)
        with open(self.trip_out_json2, 'w') as f:
            json.dump(species_triplet2, f, indent=2)

        log(f"Saved mutation counts to {self.out_json1} and {self.out_json2}", self.verbose)
        log(f"Saved triplet counts to {self.trip_out_json1} and {self.trip_out_json2}", self.verbose)


    @staticmethod
    def parse_line(line):
        parts = line.strip().split('\t')
        return parts[:5] + parts[6:-1] if len(parts) >= 9 else None

    @staticmethod
    def get_nuc(nuc_field):
        cleaned = nuc_field.translate(REMOVE_CHARS)
        return cleaned[0].upper() if cleaned else 'N'

    def extract_triplets(self, fields_list):
        context = [[], [], []]
        for fields in fields_list:
            ref_nuc = self.get_nuc(fields[REF_NUC_IDX])
            nuc1 = self.get_nuc(fields[NUC_1_IDX])
            nuc2 = self.get_nuc(fields[NUC_2_IDX])
            context[REF_IDX].append(ref_nuc)
            context[TAXA1_IDX].append(nuc1 if nuc1 not in {',', '.'} else ref_nuc)
            context[TAXA2_IDX].append(nuc2 if nuc2 not in {',', '.'} else ref_nuc)
        return context

    def detect_mutation_triplet(self, triplets):
        t1_mut = t2_mut = t1_3mer = t2_3mer = None

        if triplets[REF_IDX][PREV_IDX] == triplets[TAXA1_IDX][PREV_IDX] == triplets[TAXA2_IDX][PREV_IDX] and \
           triplets[REF_IDX][NEXT_IDX] == triplets[TAXA1_IDX][NEXT_IDX] == triplets[TAXA2_IDX][NEXT_IDX]:

            ref_base = triplets[REF_IDX][CUR_IDX]
            context = triplets[REF_IDX][PREV_IDX] + ref_base + triplets[REF_IDX][NEXT_IDX]
            if ref_base == triplets[TAXA1_IDX][CUR_IDX] and ref_base != triplets[TAXA2_IDX][CUR_IDX]:
                t2_mut = f"{context[0]}[{ref_base}>{triplets[TAXA2_IDX][CUR_IDX]}]{context[2]}"
                t1_3mer, t2_3mer = context, context
            elif ref_base == triplets[TAXA2_IDX][CUR_IDX] and ref_base != triplets[TAXA1_IDX][CUR_IDX]:
                t1_mut = f"{context[0]}[{ref_base}>{triplets[TAXA1_IDX][CUR_IDX]}]{context[2]}"
                t1_3mer, t2_3mer = context, context
            elif ref_base == triplets[TAXA2_IDX][CUR_IDX] == triplets[TAXA1_IDX][CUR_IDX]:
                t1_3mer, t2_3mer = context, context
            
        return t1_mut, t2_mut, t1_3mer, t2_3mer

    '''
    @staticmethod
    def quality_check(fields):
        return fields and '*' not in fields[NUC_1_IDX] and '*' not in fields[NUC_2_IDX] and \
               MutationExtractor.all_same(fields[NUC_1_IDX].translate(REMOVE_CHARS)) and \
               MutationExtractor.all_same(fields[NUC_2_IDX].translate(REMOVE_CHARS))
    '''
    @staticmethod
    def quality_check(fields):
        if not fields:
            return False
        nuc1 = fields[NUC_1_IDX].translate(REMOVE_CHARS).replace(',', '.').lower()
        nuc2 = fields[NUC_2_IDX].translate(REMOVE_CHARS).replace(',', '.').lower()
        return (
            '*' not in fields[NUC_1_IDX] and
            '*' not in fields[NUC_2_IDX] and
            int(fields[N_READS_1_IDX]) >= 3 and
            int(fields[N_READS_2_IDX]) >= 3 and
            MutationExtractor.all_same(nuc1) and
            MutationExtractor.all_same(nuc2)
        )

    @staticmethod
    def all_same(seq):
        return len(seq) > 0 and all(ch == seq[0] for ch in seq)


FLANK = 2

class FiveMerExtractor:
    def __init__(self, reference, taxon1, taxon2, pileup_file, output_dir, no_cache=False, verbose=True):
        self.reference = reference
        self.taxon1 = taxon1
        self.taxon2 = taxon2
        self.pileup_file = pileup_file
        self.json1_path = os.path.join(output_dir, f"{taxon1}__{taxon2}__{reference}__5mers.json")
        self.json2_path = os.path.join(output_dir, f"{taxon2}__{taxon1}__{reference}__5mers.json")
        self.no_cache = no_cache
        self.verbose = verbose
        os.makedirs(output_dir, exist_ok=True)

    def parse_line(self, line):
        parts = line.strip().split('\t')
        return parts[:5] + parts[6:-1] if len(parts) >= 9 else None

    def get_nuc(self, field):
        cleaned = field.translate(REMOVE_CHARS)
        return cleaned[0].upper() if cleaned else 'N'
    '''
    def quality_check(self, fields):
        return fields and '*' not in fields[NUC_1_IDX] and '*' not in fields[NUC_2_IDX] and \
               self.all_same(fields[NUC_1_IDX].translate(REMOVE_CHARS)) and \
               self.all_same(fields[NUC_2_IDX].translate(REMOVE_CHARS))
    '''

    def quality_check(self, fields):
        if not fields:
            return False
        nuc1 = fields[NUC_1_IDX].translate(REMOVE_CHARS).replace(',', '.').lower()
        nuc2 = fields[NUC_2_IDX].translate(REMOVE_CHARS).replace(',', '.').lower()
        return (
            '*' not in fields[NUC_1_IDX] and
            '*' not in fields[NUC_2_IDX] and
            int(fields[N_READS_1_IDX]) >= 3 and
            int(fields[N_READS_2_IDX]) >= 3 and
            self.all_same(nuc1) and
            self.all_same(nuc2)
        )


    def all_same(self, seq):
        return len(seq) > 0 and all(b == seq[0] for b in seq)

    def extract_5mer(self, fields_list):
        sequences = [[], [], []]
        for fields in fields_list:
            ref_nuc = self.get_nuc(fields[REF_NUC_IDX])
            nuc1 = self.get_nuc(fields[NUC_1_IDX])
            nuc2 = self.get_nuc(fields[NUC_2_IDX])
            nuc1 = ref_nuc if nuc1 in {',', '.'} else nuc1
            nuc2 = ref_nuc if nuc2 in {',', '.'} else nuc2
            sequences[REF_IDX].append(ref_nuc)
            sequences[TAXA1_IDX].append(nuc1)
            sequences[TAXA2_IDX].append(nuc2)
        return sequences

    def detect_mutation_5mer(self, five_mers):
        flank_indices = [0, 1, 3, 4]
        center = 2
        for i in flank_indices:
            if not (five_mers[REF_IDX][i] == five_mers[TAXA1_IDX][i] == five_mers[TAXA2_IDX][i]):
                return None, None
        ref_base = five_mers[REF_IDX][center]
        t1_base = five_mers[TAXA1_IDX][center]
        t2_base = five_mers[TAXA2_IDX][center]
        context = ''.join(five_mers[REF_IDX])
        t1_mut = f"{context[:2]}[{ref_base}>{t1_base}]{context[3:]}" if t1_base != ref_base and t2_base == ref_base else None
        t2_mut = f"{context[:2]}[{ref_base}>{t2_base}]{context[3:]}" if t2_base != ref_base and t1_base == ref_base else None
        return t1_mut, t2_mut

    def extract(self):
        if all(os.path.exists(p) for p in [self.json1_path, self.json2_path]) and not self.no_cache:
            if self.verbose:
                log("5-mer mutation files exist. Skipping.", self.verbose)
            return self.json1_path, self.json2_path

        species_mut1 = defaultdict(int)
        species_mut2 = defaultdict(int)

        with gzip.open(self.pileup_file, 'rt') as f:
            window = [None] * (2 * FLANK + 1)
            qc = [False] * (2 * FLANK + 1)
            for i in range(2 * FLANK):
                window[i] = self.parse_line(f.readline())
                qc[i] = self.quality_check(window[i])

            for line in f:
                window = window[1:] + [self.parse_line(line)]
                qc = qc[1:] + [self.quality_check(window[-1])]

                if all(qc):
                    five_mers = self.extract_5mer(window)
                    mut1, mut2 = self.detect_mutation_5mer(five_mers)
                    if mut1:
                        species_mut1[mut1] += 1
                    if mut2:
                        species_mut2[mut2] += 1

        with open(self.json1_path, 'w') as f:
            json.dump(species_mut1, f, indent=2)
        with open(self.json2_path, 'w') as f:
            json.dump(species_mut2, f, indent=2)

        log(f"Written: {self.json1_path}, {self.json2_path}", self.verbose)

        return self.json1_path, self.json2_path
    

class TripletExtractor:
    def __init__(self, reference, taxon1, taxon2, pileup_file, output_dir, no_cache=False, verbose=True):
        self.reference = reference
        self.taxon1 = taxon1
        self.taxon2 = taxon2
        self.pileup_file = pileup_file
        self.output_dir = output_dir
        self.no_cache = no_cache
        self.verbose = verbose

        os.makedirs(self.output_dir, exist_ok=True)

        self.pileup_file = pileup_file

        self.out_json1 = os.path.join(
            self.output_dir, f"{self.taxon1}__{self.taxon2}__{self.reference}__triplets.json"
        )
        self.out_json2 = os.path.join(
            self.output_dir, f"{self.taxon2}__{self.taxon1}__{self.reference}__triplets.json"
        )

    def all_same(self, seq):
        return len(seq) > 0 and all(ch == seq[0] for ch in seq)

    def get_nuc(self, field):
        cleaned = field.translate(REMOVE_CHARS)
        return cleaned[0].upper() if cleaned else 'N'

    def parse_line(self, line):
        parts = line.strip().split('\t')
        return parts[:5] + parts[6:-1] if len(parts) >= 9 else None

    def passes_qc(self, fields):
        return fields is not None and \
               '*' not in fields[NUC_1_IDX] and '*' not in fields[NUC_2_IDX] and \
               self.all_same(fields[NUC_1_IDX].translate(REMOVE_CHARS)) and \
               self.all_same(fields[NUC_2_IDX].translate(REMOVE_CHARS))

    def relevant_triplet(self, triplets):
        if triplets[REF_IDX][PREV_IDX] == triplets[TAXA1_IDX][PREV_IDX] == triplets[TAXA2_IDX][PREV_IDX] and \
        triplets[REF_IDX][NEXT_IDX] == triplets[TAXA1_IDX][NEXT_IDX] == triplets[TAXA2_IDX][NEXT_IDX]:
            if triplets[REF_IDX][CUR_IDX] == triplets[TAXA1_IDX][CUR_IDX] or \
                triplets[REF_IDX][CUR_IDX] == triplets[TAXA2_IDX][CUR_IDX]:
                return True
        return False

    def extract_triplets(self, line_fields):
        sequences = [[], [], []]
        for fields in line_fields:
            ref_nuc = self.get_nuc(fields[REF_NUC_IDX])
            nuc1 = self.get_nuc(fields[NUC_1_IDX])
            nuc2 = self.get_nuc(fields[NUC_2_IDX])
            nuc1 = ref_nuc if nuc1 in {',', '.'} else nuc1
            nuc2 = ref_nuc if nuc2 in {',', '.'} else nuc2
            sequences[REF_IDX].append(ref_nuc)
            sequences[TAXA1_IDX].append(nuc1)
            sequences[TAXA2_IDX].append(nuc2)
        return sequences

    def extract(self):
        if all(os.path.exists(p) for p in [self.out_json1, self.out_json2]) and not self.no_cache:
            log("Triplet counts already exist. Skipping.", self.verbose)
            return

        triplet_dict1 = defaultdict(int)
        triplet_dict2 = defaultdict(int)

        with gzip.open(self.pileup_file, 'rt') as f:
            line_fields = [None, self.parse_line(f.readline()), self.parse_line(f.readline())]
            qc_flags = [False, self.passes_qc(line_fields[1]), self.passes_qc(line_fields[2])]

            for line in f:
                line_fields = [line_fields[CUR_IDX], line_fields[NEXT_IDX], self.parse_line(line)]
                qc_flags = [qc_flags[CUR_IDX], qc_flags[NEXT_IDX], self.passes_qc(line_fields[NEXT_IDX])]

                if all(qc_flags):
                    triplets = self.extract_triplets(line_fields)
                    # triplet_dict1[''.join(triplets[TAXA1_IDX])] += 1
                    # triplet_dict2[''.join(triplets[TAXA2_IDX])] += 1
                    if self.relevant_triplet(triplets):
                        triplet = ''.join(triplets[REF_IDX])
                        triplet_dict1[triplet] += 1
                        triplet_dict2[triplet] += 1

        with open(self.out_json1, 'w') as f:
            json.dump(triplet_dict1, f, indent=2)
        with open(self.out_json2, 'w') as f:
            json.dump(triplet_dict2, f, indent=2)

        log(f"Triplet dictionaries written to:\n  • {self.out_json1}\n  • {self.out_json2}", self.verbose)


import os
import json
import re
from collections import defaultdict
import pandas as pd

class MutationNormalizer:
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    valid_bases = {'A', 'C', 'G', 'T'}
    mutation_pattern = re.compile(r"^[ACGT]\[[ACGT]>[ACGT]\][ACGT]$")

    def __init__(self, input_dir, output_dir=None, divergence_time=None, verbose=True):
        self.input_dir = input_dir
        self.mutation_dir = os.path.join(input_dir, "Mutations")
        self.triplet_dir = os.path.join(input_dir, "Triplets")
        self.output_dir = output_dir or os.path.join(input_dir, "Tables")
        self.divergence_time = divergence_time
        self.verbose = verbose
        os.makedirs(self.output_dir, exist_ok=True)

        self.all_collapsed_mut = {}
        self.all_norm_mut = {}
        self.all_scaled_mut = {}
        self.all_triplets = {}

    def load_json(self, path):
        with open(path, 'r') as f:
            return json.load(f)

    def save_json(self, obj, path):
        with open(path, 'w') as f:
            json.dump(obj, f, indent=2)

    def collapse_triplets(self, triplet_dict):
        collapsed = defaultdict(int)
        for triplet, count in triplet_dict.items():
            if triplet[1] in {'G', 'A'}:
                rc = ''.join(self.complement.get(nuc, nuc) for nuc in reversed(triplet))
                collapsed[rc] += count
            else:
                collapsed[triplet] += count
        return collapsed

    def get_complement(self, mutation):
        comp = [self.complement[nuc] if nuc in self.complement else nuc for nuc in mutation]
        comp[0], comp[-1] = comp[-1], comp[0]
        return ''.join(comp)

    def collapse_mutations(self, mutation_dict):
        collapsed = defaultdict(int)
        for mutation, count in mutation_dict.items():
            if mutation[2] in {'A', 'G'}:
                collapsed[self.get_complement(mutation)] += int(count)
            else:
                collapsed[mutation] += int(count)
        return collapsed

    def filter_mutations_dict(self, d):
        return {k: v for k, v in d.items() if self.mutation_pattern.match(k)}

    def filter_triplets_dict(self, d):
        return {k: v for k, v in d.items() if 'N' not in k and all(x in self.valid_bases for x in k)}

    def normalize_by_triplets(self, mutations, triplets):
        return {
            k: v / triplets.get(f"{k[0]}{k[2]}{k[-1]}", 1) if triplets.get(f"{k[0]}{k[2]}{k[-1]}", 0) > 0 else 0
            for k, v in mutations.items()
        }

    def scale_counts(self, d, target_sum=10000):
        total = sum(d.values())
        if total == 0:
            return {k: 0 for k in d}
        return {k: round(v / total * target_sum) for k, v in d.items()}

    def normalize(self):
        for file in os.listdir(self.mutation_dir):
            if not file.endswith(".json"):
                continue

            mutation_path = os.path.join(self.mutation_dir, file)
            triplet_path = os.path.join(self.triplet_dir, file.replace("mutations", "triplets"))
            key = file.replace('.json', '')

            if not os.path.exists(triplet_path):
                continue

            mutations = self.filter_mutations_dict(self.load_json(mutation_path))
            triplets = self.filter_triplets_dict(self.load_json(triplet_path))

            collapsed_mut = self.collapse_mutations(mutations)
            collapsed_tri = self.collapse_triplets(triplets)

            self.all_collapsed_mut[key] = collapsed_mut
            self.all_triplets[key] = collapsed_tri

            norm_mut = self.normalize_by_triplets(collapsed_mut, collapsed_tri)
            scaled_mut = self.scale_counts(collapsed_mut)
            scaled_norm_mut = self.scale_counts(norm_mut)

            self.all_norm_mut[key] = scaled_norm_mut
            self.all_scaled_mut[key] = scaled_mut

        self._export()

    def _export(self):
        collapsed_mutations_df = pd.DataFrame(self.all_collapsed_mut)
        collapsed_mutations_df.to_csv(os.path.join(self.output_dir, "collapsed_mutations.tsv"), sep='\t')
        pd.DataFrame(self.all_norm_mut).to_csv(os.path.join(self.output_dir, "normalized_scaled.tsv"), sep='\t')
        pd.DataFrame(self.all_scaled_mut).to_csv(os.path.join(self.output_dir, "scaled_raw.tsv"), sep='\t')
        triplets_df = pd.DataFrame(self.all_triplets)
        triplets_df.to_csv(os.path.join(self.output_dir, "triplets.tsv"), sep='\t')

        mutations_per_triplet = collapsed_mutations_df.sum() / triplets_df.sum()
        log("Mutations per triplet from divergence:", self.verbose)
        for col in mutations_per_triplet.index:
            log(f"{col.split('__')[0]}: {mutations_per_triplet[col]:.2e}", self.verbose)

        if self.divergence_time:
            mutation_rates = mutations_per_triplet / (float(self.divergence_time) * 1_000_000)
            log("Estimated mutation rates per site per year:", self.verbose)
            for col in mutation_rates.index:
                log(f"{col.split('__')[0]}: {mutation_rates[col]:.2e}", self.verbose)

