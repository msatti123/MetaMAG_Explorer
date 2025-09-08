#!/usr/bin/env python3
"""
Simple MAG Merger that properly handles existing taxonomy

This script:
1. Loads your existing nodes.dmp and names.dmp
2. Understands what taxa already exist
3. Adds new MAGs without creating duplicates
4. Handles conflicts consistently
"""

import os
import hashlib
from collections import defaultdict

def parse_existing_taxonomy(nodes_file, names_file):
    """Load and understand existing taxonomy"""
    
    print("Loading existing taxonomy...")
    
    # Parse nodes.dmp
    nodes = {}  # taxid -> {parent, rank}
    with open(nodes_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t|\t')
            taxid = int(parts[0].strip())
            parent = int(parts[1].strip())
            rank = parts[2].strip()
            nodes[taxid] = {'parent': parent, 'rank': rank}
    
    # Parse names.dmp (scientific names only)
    names = {}  # taxid -> name
    name_to_taxid = {}  # name -> taxid
    with open(names_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t|\t')
            taxid = int(parts[0].strip())
            name = parts[1].strip()
            name_type = parts[3].strip().rstrip('\t|')
            
            if name_type == 'scientific name':
                names[taxid] = name
                name_to_taxid[name] = taxid
    
    # Build taxonomy paths (taxid -> full lineage)
    print("Building taxonomy paths...")
    paths = {}
    
    def get_path(taxid):
        if taxid in paths:
            return paths[taxid]
        
        path = []
        current = taxid
        while current != 1 and current in nodes:
            if current in names:
                path.append(names[current])
            current = nodes[current]['parent']
        
        path.reverse()
        paths[taxid] = path
        return path
    
    # Build all paths
    for taxid in nodes:
        get_path(taxid)
    
    # Find max taxid
    max_taxid = max(nodes.keys())
    
    print(f"Loaded {len(nodes)} taxa, max taxID: {max_taxid}")
    
    return nodes, names, name_to_taxid, paths, max_taxid

def find_existing_taxon(taxon_name, parent_taxid, nodes, names, name_to_taxid):
    """Check if this exact taxon already exists"""
    
    if taxon_name in name_to_taxid:
        existing_taxid = name_to_taxid[taxon_name]
        # Check if parent matches
        if nodes[existing_taxid]['parent'] == parent_taxid:
            return existing_taxid
    
    return None

def generate_unique_name(taxon_name, parent_name):
    """Generate deterministic unique name for conflicts"""
    
    # Create stable hash
    key = f"{taxon_name}|{parent_name}"
    hash_suffix = hashlib.md5(key.encode()).hexdigest()[:4]
    
    # Shorten parent name for readability
    parent_short = parent_name.replace('_', '')[:4]
    
    return f"{taxon_name}_{parent_short}{hash_suffix}"

def process_new_mags(taxonomy_file, nodes, names, name_to_taxid, paths, start_taxid):
    """Process new MAGs and add to existing taxonomy"""
    
    print(f"\nProcessing new MAGs from {taxonomy_file}")
    
    # Track new additions
    new_nodes = {}
    new_names = {}
    new_genomes = {}
    conflicts_resolved = []
    
    # Current new taxID
    next_taxid = start_taxid
    
    ranks = ('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species')
    
    with open(taxonomy_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split('\t')
            if len(parts) != 2:
                continue
            
            genome_id, lineage = parts
            
            # Skip header
            if genome_id.lower() in ['genome_id', 'user_genome']:
                continue
            
            # Parse lineage
            taxa = lineage.split(';')
            if len(taxa) != 7:
                print(f"  Warning: Skipping incomplete lineage at line {line_num}")
                continue
            
            # Process each level
            parent_taxid = 1
            parent_name = 'root'
            
            for i, taxon in enumerate(taxa):
                if '__' not in taxon:
                    break
                
                code, taxon_name = taxon.split('__', 1)
                if not taxon_name:
                    break
                
                rank = ranks[i]
                
                # Check if exists with this exact parent
                existing = find_existing_taxon(taxon_name, parent_taxid, nodes, names, name_to_taxid)
                
                if existing:
                    # Reuse existing
                    current_taxid = existing
                    current_name = taxon_name
                else:
                    # Check if name exists under different parent
                    if taxon_name in name_to_taxid:
                        # Conflict! Generate unique name
                        unique_name = generate_unique_name(taxon_name, parent_name)
                        
                        # Check if this unique name already exists
                        if unique_name in name_to_taxid:
                            current_taxid = name_to_taxid[unique_name]
                            current_name = unique_name
                        else:
                            # Create new with unique name
                            current_taxid = next_taxid
                            next_taxid += 1
                            
                            new_nodes[current_taxid] = {'parent': parent_taxid, 'rank': rank}
                            new_names[current_taxid] = unique_name
                            
                            # Update lookups
                            nodes[current_taxid] = new_nodes[current_taxid]
                            names[current_taxid] = unique_name
                            name_to_taxid[unique_name] = current_taxid
                            
                            conflicts_resolved.append({
                                'original': taxon_name,
                                'renamed': unique_name,
                                'rank': rank
                            })
                            
                            current_name = unique_name
                    else:
                        # Truly new taxon
                        current_taxid = next_taxid
                        next_taxid += 1
                        
                        new_nodes[current_taxid] = {'parent': parent_taxid, 'rank': rank}
                        new_names[current_taxid] = taxon_name
                        
                        # Update lookups
                        nodes[current_taxid] = new_nodes[current_taxid]
                        names[current_taxid] = taxon_name
                        name_to_taxid[taxon_name] = current_taxid
                        
                        current_name = taxon_name
                
                # Update for next iteration
                parent_taxid = current_taxid
                parent_name = current_name
            
            # Store genome mapping
            new_genomes[genome_id] = parent_taxid
    
    print(f"  Added {len(new_nodes)} new taxa")
    print(f"  Added {len(new_genomes)} new genomes")
    print(f"  Resolved {len(conflicts_resolved)} conflicts")
    
    return new_nodes, new_names, new_genomes, conflicts_resolved

def write_merged_taxonomy(all_nodes, all_names, new_genomes, conflicts, output_dir):
    """Write the merged taxonomy files"""
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Write nodes.dmp
    print("\nWriting merged nodes.dmp...")
    with open(os.path.join(output_dir, 'nodes.dmp'), 'w') as f:
        for taxid in sorted(all_nodes.keys()):
            info = all_nodes[taxid]
            f.write(f"{taxid}\t|\t{info['parent']}\t|\t{info['rank']}\t|\t\t|\n")
    
    # Write names.dmp
    print("Writing merged names.dmp...")
    with open(os.path.join(output_dir, 'names.dmp'), 'w') as f:
        for taxid in sorted(all_names.keys()):
            name = all_names[taxid]
            f.write(f"{taxid}\t|\t{name}\t|\t\t|\tscientific name\t|\n")
    
    # Write taxid.map for new genomes
    print("Writing taxid.map for new genomes...")
    with open(os.path.join(output_dir, 'taxid.map'), 'w') as f:
        for genome_id, taxid in sorted(new_genomes.items()):
            f.write(f"{genome_id}\t{taxid}\n")
    
    # Write conflicts log
    if conflicts:
        print("Writing conflicts.log...")
        with open(os.path.join(output_dir, 'conflicts.log'), 'w') as f:
            f.write("Resolved Naming Conflicts\n")
            f.write("=" * 50 + "\n\n")
            for c in conflicts:
                f.write(f"{c['rank']}: {c['original']} -> {c['renamed']}\n")

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Merge new MAGs with existing taxonomy')
    parser.add_argument('--existing-nodes', required=True, help='Existing nodes.dmp')
    parser.add_argument('--existing-names', required=True, help='Existing names.dmp')
    parser.add_argument('--new-taxonomy', required=True, help='New MAGs taxonomy.tsv')
    parser.add_argument('--output-dir', default='merged_taxonomy', help='Output directory')
    
    args = parser.parse_args()
    
    # Load existing taxonomy
    nodes, names, name_to_taxid, paths, max_taxid = parse_existing_taxonomy(
        args.existing_nodes, args.existing_names
    )
    
    # Process new MAGs
    new_nodes, new_names, new_genomes, conflicts = process_new_mags(
        args.new_taxonomy, nodes, names, name_to_taxid, paths, max_taxid + 1
    )
    
    # Write merged taxonomy
    write_merged_taxonomy(nodes, names, new_genomes, conflicts, args.output_dir)
    
    print(f"\nDone! Merged taxonomy in {args.output_dir}/")
    print(f"Total taxa: {len(nodes)}")
    print(f"Total genomes added: {len(new_genomes)}")

if __name__ == '__main__':
    main()