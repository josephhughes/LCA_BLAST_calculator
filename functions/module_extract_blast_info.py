#! /usr/bin/env python3
import collections
from tqdm import tqdm
import os

def format_info(format):
    formatting = format.split('+')
    count = 0
    pident = 'NA'
    qcovs = 'NA'
    staxid = 'NA'
    qaccver = 'NA'
    for item in formatting:
        if item == 'pident':
            pident = count
        elif item == 'qcovs':
            qcovs = count
        elif item == 'staxid':
            staxid = count
        elif item == 'qaccver':
            qaccver = count
        count = count + 1

    return pident, qcovs, staxid, qaccver
        
def extract_blast(input, pident, qcovs, staxid, qaccver):
    LCA_dict = collections.defaultdict(list)
    taxid_list = []
    with tqdm(total = os.path.getsize(input)) as pbar:
        with open(input, 'r') as infile:
            for line in infile:
                pbar.update(len(line))
                cols = line.strip().split('\t')
                if len(cols) <= max(pident, qcovs, staxid, qaccver):
                    continue
                pident_val = cols[pident].strip()
                qcovs_val = cols[qcovs].strip()
                staxid_val = cols[staxid].strip()
                qaccver_val = cols[qaccver].strip()
                LCA_dict[qaccver_val].append([staxid_val, pident_val, qcovs_val])
                taxid_list.append(staxid_val)
    taxid_set = set(taxid_list)

    return LCA_dict, taxid_set

def extract_nodes(nodes):
    taxids = {}
    with tqdm(total = os.path.getsize(nodes)) as pbar:
        with open(nodes, 'r') as infile:
            for line in infile:
                pbar.update(len(line))
                tax = line.split('\t|\t')[0]
                taxup = line.split('\t|\t')[1]
                rank = line.split('\t|\t')[2]
                taxids[tax] = [rank, taxup]

    return taxids


def extract_names(names):
    name_dict = {}
    with tqdm(total = os.path.getsize(names)) as pbar:
        with open(names) as infile:
            for line in infile:
                pbar.update(len(line))
                sc = line.split('\t')[6]
                if sc == 'scientific name':
                    taxid_name = line.split('\t|\t')[0]
                    name = line.split('\t|\t')[1].replace(' ', '_')
                    name_dict[taxid_name] = name

    return name_dict

def tax_lineage(taxid_set, ranks, taxid_info, names_info):
    """
    Generate taxonomic lineages for a set of TaxIDs.
    Uses caching to avoid redundant tree traversals.
    """
    rank_list = ranks.split('+')
    rank_dict = {rank: True for rank in rank_list}
    true_lineage = collections.defaultdict(list)
    
    # Cache for full lineage paths to avoid re-traversing the same nodes
    # Stores: {taxid: {rank: taxid, ...}}
    path_cache = {}

    for tax in tqdm(taxid_set, desc="Resolving lineages"):
        ktax = tax
        if ktax not in taxid_info:
            true_lineage[ktax] = [[rank, 'nan', 'nan'] for rank in rank_list]
            continue
        
        # Build full taxonomic path for this TaxID
        tax_path = {}
        curr = tax
        visited = set()
        
        while curr in taxid_info:
            if curr in visited: # Cycle detection
                break
            visited.add(curr)
            
            rank_found = taxid_info[curr][0]
            # Map 'domain' to 'superkingdom' if needed
            if rank_found == 'domain' and 'superkingdom' in rank_dict:
                rank_found = 'superkingdom'
            
            if rank_found in rank_dict:
                tax_path[rank_found] = curr
            
            parent = taxid_info[curr][1]
            if curr == parent: # Reached root
                break
            
            # Use cache if we've already resolved the rest of the path
            if parent in path_cache:
                tax_path.update(path_cache[parent])
                break
                
            curr = parent
            
        # Store resolved path in cache (for future use by other TaxIDs)
        path_cache[ktax] = tax_path
        
        # Format lineage according to requested ranks
        for rank in rank_list:
            v = tax_path.get(rank, 'nan')
            name = names_info.get(v, 'nan') if v != 'nan' else 'nan'
            true_lineage[ktax].append([rank, v, name])
        
    return true_lineage

def calculate_lca(LCA_dict, PIDENT, QCOV, RANKS, lineage):
    """
    Calculate the Lowest Common Ancestor for each query based on filtered BLAST hits.
    Uses list-based prefix matching for robustness.
    """
    LCA_final = {}
    ranks_list = RANKS.split('+')
    
    for item in tqdm(LCA_dict, desc="Calculating LCA"):
        hits_lineages = []
        blast_hit_number = 0
        
        for v in LCA_dict[item]:
            # Filter by percent identity and query coverage
            if float(v[1]) >= float(PIDENT) and float(v[2]) >= float(QCOV):
                if v[0] not in lineage:
                    continue
                
                # Get the lineage as a list of scientific names
                names = [rank_info[2] for rank_info in lineage[v[0]]]
                hits_lineages.append(names)
                blast_hit_number += 1
        
        if not hits_lineages:
            LCA_final[item] = {
                'lineage': 'na',
                'rank': 'na',
                'hits': 0
            }
            continue

        # Find common prefix across all hits
        common_names = []
        for i in range(len(ranks_list)):
            level_names = set(h[i] for h in hits_lineages)
            if len(level_names) == 1 and 'nan' not in level_names:
                common_names.append(list(level_names)[0])
            else:
                break
        
        result_lineage = ";".join(common_names) if common_names else 'NA'
        obtained_rank = ranks_list[len(common_names)-1] if common_names else 'NA'
        
        LCA_final[item] = {
            'lineage': result_lineage,
            'rank': obtained_rank,
            'hits': blast_hit_number
        }
    
    return LCA_final
                
def read_otutable(FREQ):
    otutable = {}
    with tqdm(total = os.path.getsize(FREQ)) as pbar:
        with open(FREQ, 'r') as infile:
            for line in infile:
                pbar.update(len(line))
                line = line.rstrip('\n')
                otu_name = line.split('\t', 1)[0]
                samples = line.split('\t', 1)[1]
                otutable[otu_name] = samples
    
    return otutable

def combine_data(otutable, LCA_final, OUTPUT):
    """
    Combine LCA results with the original OTU frequency table.
    """
    print(f'writing combined data to {OUTPUT}')
    with open(OUTPUT, 'w') as outfile:
        # Sort and process in order
        for item in tqdm(otutable, desc="Exporting"):
            if item.startswith('#'):
                # Header row: add our new columns
                headers = [
                    'qaccver',
                    'rank',
                    'blast_hit_number',
                    'taxonomic_lineage',
                    otutable[item]
                ]
                outfile.write("\t".join(headers) + "\n")
                continue
            
            samples = otutable[item]
            if item not in LCA_final:
                lineage = 'na'
                rank = 'na'
                hits = 'na'
            else:
                res = LCA_final[item]
                lineage = res['lineage']
                rank = res['rank']
                hits = str(res['hits'])
            
            row = [item, rank, hits, lineage, samples]
            outfile.write("\t".join(row) + "\n")

def export_lca(LCA_final, output):
    """
    Export LCA results to a tab-delimited file.
    """
    print(f'exporting LCA results to {output}')
    with open(output, 'w') as outfile:
        # Write header
        outfile.write('#qaccver\trank\tblast_hit_number\ttaxonomic_lineage\n')
        
        for item in tqdm(LCA_final, desc="Exporting"):
            res = LCA_final[item]
            row = [
                item,
                res['rank'],
                str(res['hits']),
                res['lineage']
            ]
            outfile.write("\t".join(row) + "\n")