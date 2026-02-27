import os
import csv
import pandas as pd
import numpy as np
from collections import defaultdict
import ast


def load_gc_ratios(gc_path):
    gc_map = {}
    with open(gc_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            gc_map[row['tag']] = (round(float(row['gc_ratio']),2)) * 100  # Convert to %
    return gc_map

def load_labels(label_path):
    label_map = {}
    with open(label_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            label_map[row['tag']] = int(row['label'])
    return label_map

def load_conservation(conservation_path):
    cons_map = {}
    with open(conservation_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            key = (row['gene1_id'], row['gene2_id'])
            cons_map[key] = round(float(row['conservation_score']),2)
    return cons_map


def codon_bias_distance(peg1, peg2, codon_bias_lookup):
    """
    Compute the Euclidean distance between codon bias vectors for peg1 and peg2.
    Returns a float (or None if one of the pegs is missing).
    """
    vec1 = codon_bias_lookup.get(peg1)
    vec2 = codon_bias_lookup.get(peg2)
    
    if vec1 is None or vec2 is None:
        return None

    # Ensure both have all 64 codons (fill missing with 0.0)
    codons = sorted(set(vec1.keys()) | set(vec2.keys()))
    v1 = np.array([vec1.get(codon, 0.0) for codon in codons])
    v2 = np.array([vec2.get(codon, 0.0) for codon in codons])

    return np.linalg.norm(v1 - v2)


def parse_genome(path):
   
        # Read while handling null bytes and preserving column names
        df = pd.read_csv(path, sep="\t", encoding="utf-8")
        

        genes = []
        for index, row in df.iterrows():#reader:
            if row['feature_type'] != 'CDS':
                continue
            start = int(row['start'])
            peg = row['patric_id']
            length = abs(int(row['end']) - int(row['start'])) + 1
            strand = row['strand']
            product = row['product'] if row['product'] else 'hypothetical protein'
            family = row['pgfam_id'] if row['pgfam_id'] else ''
            genes.append((start, peg[4:], length, strand, product, family))
        genes.sort(key=lambda x: x[0])
        print(len(genes))
        return genes



def load_string_scores(path):
    """
    Load STRING scores from a CSV into a dictionary mapping (pegA, pegB) → score.
    Adds both (A, B) and (B, A) for symmetry.
    """
    df = pd.read_csv(path)  # This assumes header = ["pegA", "pegB", "score"]
    
    string_scores = {}
    for _, row in df.iterrows():
        pegA = row["peg1"]
        pegB = row["peg2"]
        score = int(row["string"])
        
        # Store both directions
        string_scores[(pegA, pegB)] = score
        string_scores[(pegB, pegA)] = score

    return string_scores

def create_textual_dataset(genome_dir, gc_path, label_path, conservation_path, output_path, codon = None, sscore=None):
    gc_map = load_gc_ratios(gc_path)
    label_map = load_labels(label_path)
    cons_map = load_conservation(conservation_path)


    if codon:

        codon_bias_df = pd.read_csv(codon)
        codon_bias_lookup = {
            row["tag"]: ast.literal_eval(row["bias"])
            for _, row in codon_bias_df.iterrows()
        }

    strings = None
    if sscore:
        strings = load_string_scores(sscore)
    rows = []
    f = 0
    nf = 0
    for filename in os.listdir(genome_dir):
        if not filename.endswith('.features.tab'):
            continue

        genome_id = filename.split('.')[0]
        path = os.path.join(genome_dir, filename)
        genes = parse_genome(path)

        for i in range(len(genes) - 1):
            g1 = genes[i]
            g2 = genes[i+1]
            peg1, peg2 = g1[1], g2[1]


            if peg1 not in label_map:              
                continue


            gc1 = gc_map[peg1]
            gc2 = gc_map[peg2]
            gc_diff = round(abs(gc1 - gc2), 1)

            if label_map[peg1] == 1:
                label = 1 
            elif label_map[peg1] == 0:
                label = 0
            cons_score = cons_map.get((g1[5], g2[5]), 0.0)
            if cons_score == 0:
                f += 1
            else:
                nf += 1

            intergenic = genes[i+1][0] - (genes[i][0] + g1[2])
            strand_relation = "same strand" if g1[3] == g2[3] else "different strands"

            text = (
                f"Gene A is {g1[2]} bp long, annotated as \"{g1[4]}\", and is part of the {g1[5]} family. "
                f"It is followed by Gene B, {g2[2]} bp long, annotated as \"{g2[4]}\", and is part of the {g2[5]} family. "
                f"The intergenic distance is {intergenic} base pairs. "
                f"The GC content difference is {gc_diff}%. "
                f"They are on the {strand_relation}. "
                f"They are adjacent in {cons_score*100:.1f}% of representative genomes."
                f"Their String score is unavailable."
            )

            if strings:
                if (peg1,peg2) in strings:
                    s = f"Their String score is {strings[(peg1,peg2)]}."
                    f += 1
                else:
                    s = f"Their String score is unavailable."
                    nf += 1
                stext = f"{text} {s}"
            else:
                stext = text


            # Load and parse codon bias file
            if codon:
                bias_diff = codon_bias_distance(peg1, peg2, codon_bias_lookup)
                if bias_diff is not None:
                    codon_text = f"The codon bias difference between the genes is {bias_diff:.2f}."
                else:
                    codon_text = "Codon bias information is not available for both genes."

                # Append to description
                final_description = f"{stext} {codon_text}"
                stext = final_description
            else:
                final_description = stext

            rows.append((genome_id, peg1, peg2, stext, label))

    with open(output_path, 'w', newline='') as out:
        writer = csv.writer(out)
        writer.writerow(['genome_id', 'gene1', 'gene2', 'text', 'label'])
        writer.writerows(rows)

    print(f"Wrote {len(rows)} gene pairs to {output_path}")
    print(f, nf)

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 6:
        print("Usage: python3 serialize_gene_pairs.py <genome_dir> <gc.csv> <labels.csv> <conservation.tsv> <codon?> <codon+string?> <output.csv>")
        exit(1)
    if len(sys.argv) == 8:
        create_textual_dataset(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[7], sys.argv[5], sys.argv[6])
    elif len(sys.argv) == 7:
        create_textual_dataset(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[6], sys.argv[5])
    else:
        create_textual_dataset(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
