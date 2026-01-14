import pyham
import os 
import argparse
import gffutils
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import logging 
import time 
from pyfaidx import Fasta
import json
import urllib3
import numpy as np
import pandas as pd

def build_arg_parser():
    """Handle the parameter sent when executing the script from the terminal

    Returns:
    -----------
    A parser object with the chosen option and parameters"""

    parser = argparse.ArgumentParser(description="Extract consensus sequence from a set of GFF file containing different annotation sets for a single genome and an OMA run with these annotations.")
    
    subparsers = parser.add_subparsers(title="Commands")
    subparsers.required = True
    prepare_parser = subparsers.add_parser("prepare_data", 
            help="Handles GFF file and set up the OMA database.", 
            description="Uses GFF file to produce FASTA files (and splice files as needed), then copies them into the corresponding OMA Folder.")
    prepare_parser.set_defaults(func=prepare_data)
    extract_consensus_parser = subparsers.add_parser("extract_consensus", 
            help="Uses the OMA Standalone parser to create a consensus annotation.", 
            description='Takes as input the orthoXML and species tree from OMA and creates a consensus GFF and FASTA file.')
    extract_consensus_parser.set_defaults(func=extract_consensus)   
    for subp in [prepare_parser, extract_consensus_parser]:
            subp.add_argument('-a', '--gff_annot_folder', required=True, help="Path to the folder containing the source GFF annotations.")
            subp.add_argument('-f', '--fasta_folder', required=True, help="Path to the output folder that will contain FASTA sequences corresponding to each annotation." )
    prepare_parser.add_argument('-d', '--OMA_folder', required=True, help="Path to the OMA Standalone Folder downloaded through OMA Browser Export AllvsAll.")
    prepare_parser.add_argument('-s', '--splice_folder', required=True, help="Path to the output folder that will contain splice files corresponding to each annotation in which splice variants are explicitely included.")
    prepare_parser.add_argument('-g', '--genome_file', required=True, help="Path to the genome file to be annotated, in the FASTA format.")
    prepare_parser.add_argument('-t', '--feature_type', help="Path to a tsv file indicating which feature should be used in GFF file to define gene, transcript or CDS. 4 columns by line in order: filename, gene, transcript, CDS.")

    extract_consensus_parser.add_argument('-x', '--orthoxml', required=True,  help='Path to the orthoXML file from the OMA Standalone run.')
    extract_consensus_parser.add_argument('-st', '--species_tree',required=True,  help='Path to the species tree used in the OMA Standalone run, in newick format. A copy can be found in the output folder of OMA as ManualSpeciesTree.nwk')
    extract_consensus_parser.add_argument('-o', '--output_prefix', required=True, help='Output file path and prefix. Two output files will be created: a GFF annotation file and a FASTA file.')
    extract_consensus_parser.add_argument('-t', '--feature_type', help='"Path to a tsv file indicating which feature should be used in GFF file to define gene, transcript or CDS. 4 columns by line in order: filename, gene, transcript, CDS.')
    extract_consensus_parser.add_argument('-p', '--priorities', default=None, help='Comma-separated string indicating, in order of priorities, which source to use to make the reference sequence for each consensus gene. Should contain all file prefixes. Ex: rna,abinitio,homology')

    
    return parser


def timer_wrapper(func):
    """A timer to evaluate the time any input function take to complete when called. Used for debugging.
    Args:
        funcr (function) : The function to time
    """
    def timer(*args, **kwargs): 
        """Implementation of the timer. Return the result of the function on which it was call but log the time elapsed."""
        t0 = time.time()
        result = func(*args, **kwargs)
        t1 = time.time()
        elapsed = t1 - t0
        logging.debug(f"@timer:{func.__name__} ran in {elapsed:0.4f} seconds")
        return result
    return timer

def extract_consensus(args):
    """Extract the consensus genes as obtained through the HOGs. This is the main fuction which all the different subfunction of the program
    Args:
        args : the command line arguments as provided by the user. Expect fasta_folder, gff_annot_folder, orthoxml, species_tree, output_prefix and feature_type to have value.
    """
    #Obtain parameters
    input_fasta_folder=args.fasta_folder
    input_gff_folder=args.gff_annot_folder
    orthoxml=args.orthoxml
    species_tree=args.species_tree
    output_prefix=args.output_prefix
    feature_type_file = args.feature_type
    priorities = args.priorities
    #If a feature_type_file was provided, obtain the mapping. Is used to read non-standard GFF definitions
    if feature_type_file:
        logging.info(f'A feature type file was provided. Reading {feature_type_file}')
  
        feature_type_map = extract_feature_map(feature_type_file)
    else:
        feature_type_map = None

    #List the files in fasta and gff directories
    fasta_files = [os.path.join(input_fasta_folder,f) for f in  os.listdir(input_fasta_folder)]
    gff_files = [os.path.join(input_gff_folder,f) for f in os.listdir(input_gff_folder)]
    if priorities:
        prefixes = ['.'.join(os.path.basename(x).split('.')[0:-1]) for x in gff_files]
        corresponding_file = 0
        for src in priorities.split(','):
            if src in prefixes:
                corresponding_file+=1
        if corresponding_file != len(gff_files):
            logging.info('List in priorities does not correspond to input file prefixes. Exiting.')
            exit()
    logging.info('Associating protein sequences to gene annotaton in the GFF files.')
    #Get the correspondance between genes/proteins in the GFF file and sequences in the FASTA file
    fasta_corr, gff_corr = map_gff_fasta(fasta_files, gff_files, feature_type_map)
    #Obtain the protein identifier for each "consensus" HOG
    logging.info('Extracting orthologous groups info from OMA.')
    protids, hoglist = get_protid_per_groups(species_tree, orthoxml, '/'.join(['.'.join(f.split('.')[0:-1]) for f in os.listdir(input_fasta_folder)]))
    #Select one sequence for each consensus HOG and obtain its coordinate
    logging.info('Selecting consensus annotation for each group.')
    cons_fasta, cons_gff, selected_by_src = select_consensus_sequence(protids,gff_corr, fasta_corr,priorities)
    #Get the id of selected genes, order should be the same as hoglist
    ordered_ids = [seq.id for seq in cons_fasta]
    report_data_counter, support_matrix = make_support_analysis(hoglist,ordered_ids)
    report_data = {"gene_nr": len(cons_fasta), 'support': report_data_counter, 'selected' : selected_by_src}
    #Write the output files
    logging.info('Writing output files.')
    write_report(output_prefix+'.report.txt', report_data)
    write_matrix(output_prefix+'.detailed_report.txt', support_matrix)
    write_fasta(output_prefix+'.fa', cons_fasta)
    write_gff(output_prefix+'.gff', cons_gff)


def get_protid_per_groups(tree_path, orthoxml_path, ancestor):
    """Obtain a list of list, where each sublist correspond to a "consensus HOG" and its content each of the protein ID whitin this group
    Args:
        tree_path (str): Path to the species tree
        orthoxml_path (str) : Path to the OMA OrthoXML file
        ancestor (str) : name of the ancestral node of all annotations used as species for OMAnnotation step 2
    Returns:
        hog_protid (list) : a list of list, one by HOG, that contains protein ids (str) 
        """
    hog_protid = list()
    hog_list = list()
    #Access the OrthoXML content with pyham
    tree = pyham.utils.get_newick_string(tree_path, type="nwk")
    ham_analysis = pyham.Ham(tree, orthoxml_path, use_internal_name=False)
    #The ancestral node name is defined as the concatenation of all children (annotations). This make sure the order is irrelevant.
    for ag in [g.name for g in ham_analysis.get_list_ancestral_genomes()]:
        if sorted(ag.split('/'))== sorted(ancestor.split('/') ):
            ancestor = ag
            break
    #Obtain the consensus ancestor and HOGs and create the list
    ancestor_genome = ham_analysis.get_ancestral_genome_by_name(ancestor)
    ancestral_genes = ancestor_genome.genes
    for gene in ancestral_genes:
        
        hog_protid.append([(g.prot_id,g.genome.name) for g in gene.get_all_descendant_genes()])
        hog_list.append(gene)
    return hog_protid, hog_list

def make_support_analysis(hoglist, select_seq_id):
    """Extract information about the sources supporting each of the selected genes in the dataset, as well as the sum.
    Args:
        hoglist (list): List of HOGs at the target node ('ancestor' of all of the input sources)
        select_seq_id (list) : List of identifiers of all selected genes
    Returns:
        counter (dict) : A dictionary made to count genes with support from each source and whether or not they have external homology support for the species in OMA.
        df_support (dataframe) : a gene by gene support matrix, mark from which source annotations and each of the OMA species supporting it.
        """
    counter = dict()
    support_occurence = dict()
    for index, hog in enumerate(hoglist):
        rhog = hog.get_top_level_hog()
        level_rhog = rhog.genome.name
        all_support = set([g.genome.name for g in rhog.get_all_descendant_genes()])
        sources = [l.name for l in hog.genome.taxon.get_leaves()]
        source_support = set([g.genome.name for g in hog.get_all_descendant_genes()])
        homology_support = all_support.difference(sources)
        if len(homology_support) != 0:
            source_support = source_support |{"other_species_support"}
        support_key = tuple(source_support)
        counter[support_key]= counter.get(support_key,0)+1
        for supp in all_support:
            support_occurence[supp]= support_occurence.get(supp, []) + [index]
    support_matrix = np.zeros((len(hoglist), len(support_occurence.keys())))
    for j, key in enumerate(support_occurence):
        ind_pos = support_occurence[key]
        for k in ind_pos:
            support_matrix[k][j]=1
    support_matrix = np.vstack((support_matrix,np.sum(support_matrix,axis=0)), )
    df_support = pd.DataFrame(support_matrix, columns=list(support_occurence.keys()), index=select_seq_id+['Total'])
    
    return counter, df_support       


def map_gff_fasta(fasta_file, gff_file, gff_map =None):
    """Obtain the correspondance between protein sequences in the FASTA file and the corresponding entry in the corresponding file (shared prefix).
    Args:
        fasta_file (list) : list of path (str) to all fasta files. Must be the same length as gff_file
        gff_file (list) : list of path (str) to all gff files. 
        gff_map (dict) : whether the GFF have special feature type to use instead of the standard ones. Use gff_file path as key and map as value. Default to None
    Returns:
        record_dict (dict) : a dictionary of all FASTA sequences using a contactenation of protein ID and source a key and record object as value.
        gff_dic (dict) : a dictionary of protein to GFF segment using a contactenation of protein ID and source a key and record object as value.
    """
    gff_dic = {}
    records = []
    record_dict = {}
    ordered_ffile = sorted(fasta_file)
    ordered_gfile = sorted(gff_file)
    #Make sure each FASTA file and GFF file have a one to one correspondance
    if len(ordered_ffile)!= len(ordered_gfile):
        raise ValueError('FASTA and GFF files do not correspond')

    for i, ffile in enumerate(ordered_ffile):
        gfile = ordered_gfile[i]
        prefix = '.'.join(os.path.basename(ffile).split('.')[0:-1])
        prefix_gff = '.'.join(os.path.basename(gfile).split('.')[0:-1])
        if prefix!=prefix_gff:
            raise ValueError('FASTA and GFF file do not share prefixes')
        #Use feature map if provided
        if gff_map:
            type_map = gff_map.get(gfile, False)
            if type_map:
                type_map = {v : k for k,v in type_map.items()}
        else:
            type_map=False
        #Extract FASTA to GFF correspondance and add to the dictionary
        gff_dic = {**gff_dic,**FASTA_gff_mapper(ffile, gfile, type_map)}
        #Create a dictionary of protein ID to sequence
        records = list(SeqIO.parse(ffile, 'fasta'))
        for record in records:
            record.description = record.description+' origin='+prefix
            record_dict["_".join([record.id,prefix])] = record
    return record_dict, gff_dic

def extract_FASTA_sequences(fasta_file):
    """
    Read a FASTA file and extract a dictionary of protein_id (str) : record using BioPython
    Args:
        fasta_file (str) : path to the FASTA file
    Returns:
        record_dict (dict) : a dictionary of protein_id (str) : record
    """
    record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    return record_dict

@timer_wrapper
def extract_gff(gff_file):
    """Read a GFF file and store it in memory using the gffutils library.
    Args:
        gff_file (str) : path to the GFF file
    Returns:
        db : a gffutils db object representing the GFF"""
    db = gffutils.create_db(gff_file, ':memory:', merge_strategy="create_unique", keep_order = True, checklines=1,force_gff=True)
    return db

def FASTA_gff_mapper(fasta_file, gff_file, type_map=None):
    """Create a mapping between each protein in a FASTA file and its corresponding lines in a GFF file.
    Args:
        fasta_file (str) : path to the FASTA file
        gff_file (str) : path to the GFF file
        type_map (dict) : a mapping of feature name used in the GFF to the standard name. Used to process non-standard GFF files. Default to None for standard files.
    Returns: 
        map_id_gff (dict) : a dictionary of (protein_id (str), prefix of the GFF file (str) : GFF segment (str). Used to associate each protein to its GFF.
    """
    gff_db = extract_gff(gff_file)
    prefix = '.'.join(os.path.basename(gff_file).split('.')[0:-1])
    map_id_gff = dict()
    records = extract_FASTA_sequences(fasta_file)
    #Query the GFF db for the sequence identifier, and obtain the element and its parents  (get the top level parents and all its children if needed} --> all gene informations
    for seq_id in records.keys():
        try:
            feat = gff_db[seq_id]
        except gffutils.FeatureNotFoundError:
            continue
        #Chose whether to order elements in ascending or descending position depending of the strand
        if feat.strand=='+':
            ordering = 'start'
            reverse =False
        else:
            ordering = 'end'
            reverse= True
        parents = list(gff_db.parents(seq_id))
        
        if len(parents)> 0:
            parent = parents[0]
            children = list(gff_db.children(parent.id, order_by=ordering, reverse=reverse))
        else:
            parent = feat
            children = list(gff_db.children(seq_id, order_by=ordering, reverse=reverse))
        complete_gff = [parent] +children
        gff_str= str()
        #Write the file from the obtained components
        for lines in complete_gff:
            if type_map:
                lines.featuretype = type_map.get(lines.type_map, lines.type_map)
            gff_str += str(lines)+'\n'
        map_id_gff["_".join([seq_id,prefix])]= gff_str
    return map_id_gff


def get_longest_seq(seq_list):
    """Select the longest sequence in a list of BioPython records.
    Args:
        seq_list(list) : a list of records
    Returns:
        out_index (int) : the index of the selected record
        out_record (Bio.Record) : the selected record
    """
    max_seq_len = 0
    for index, record in enumerate(seq_list):
        #Select the longuest sequence seen in the list, in case of tie select the first it sees
        if len(record.seq) > max_seq_len:
            max_seq_len = len(record.seq)
            out_record = record
            out_index = index
    return out_index, out_record

def get_seq_by_prio(hog_prot_list, priorities,corr_fasta_map):
    """Select the the representative sequence according to the priority of the source.
    Args:
        hog_prot_list(list) : a list of protein id and source
        priorities(str) : a comma-separated list of source in order of more to least trustable
        corr_fasta_map(dict) : a dict linking a couple of protein id and source to the corresponding sequence in the source FASTA file
    Returns:
        seq_id (str) : a sequence identifier composed of the sequence id and the prefix of its source
        sequence (Bio.Record) : the selected BioPython record
    """
    seq_id=None
    for source in priorities.split(','):
        for identifier, sfile in hog_prot_list:
            if sfile==source:
                seq_id = "_".join([identifier.split(' ')[0], sfile])
                sequence = corr_fasta_map[seq_id]
        if seq_id:
            break
    return seq_id, sequence


def select_consensus_sequence(hog_prot_id_list, corr_gff_map, corr_fasta_map, source_priorities=None):
    """Select the longest sequence as the best representative of any individual consensus genes (as representend by a HOG at the ancestor of annotations) for all 
    of the consensus genes. Returns all consensus sequences and GFF as lists.
    Args:
        hog_prot_id_list (list): A list of list, where each list correspond to a consensus HOG and contains all of the protein_ids that are involved in the consensus.
        corr_gff_map (dict) : a dictionary of (protein_id, source of the annotation) to a GFF string.
        corr_fasta_map (dict) : a dictionary of a protein_id, source of the annotation to a BioPython record, containing the sequence
    Returns:
        consensus_seq (list) : a list of all the selected sequences
        consensus_gff (list) : a list of all the selected genes's GFF segment
    """
    consensus_seq = []
    consensus_gff = []
    chosen_src_nr = {}
    for hog_prot_id in hog_prot_id_list:
        if not source_priorities:
            hog_record_list = [corr_fasta_map["_".join([hid.split(' ')[0], hif])] for hid, hif in hog_prot_id]
            chosen_seq_index, chosen_seq = get_longest_seq(hog_record_list)
            chosen_seq_id = "_".join([hog_prot_id[chosen_seq_index][0].split(' ')[0], hog_prot_id[chosen_seq_index][1]])
        else:
            chosen_seq_id, chosen_seq = get_seq_by_prio(hog_prot_id, source_priorities,corr_fasta_map)
        chosen_src = chosen_seq_id.split('_')[-1]
        chosen_src_nr[chosen_src] = chosen_src_nr.get(chosen_src, 0)+1
        corr_gff = corr_gff_map[chosen_seq_id]
        consensus_gff.append(corr_gff)
        consensus_seq.append(chosen_seq)
    selected_by_src =  chosen_src_nr
    return consensus_seq, consensus_gff, selected_by_src

def write_fasta(ofile, consensus_seq):
    """A generalist function to write a FASTA file from a list of records using BioPython
    Args: 
        ofile (str) : Path to the output file
        consensus_seq (list) : A list of records to be  written as a FASTA file
    """
    with open(ofile, 'w') as output_handle:
        SeqIO.write(consensus_seq, output_handle, 'fasta')


def write_gff(ofile, consensus_gff):
    """A generalist function to write a GFF file from a list of GFF segments by concatenating them.
    Args: 
        ofile (str) : Path to the output file
        consensus_gff (list) : a list of GFF segments to be written together
    """
    with open(ofile, 'w') as output_handle:
        output_handle.write("".join(consensus_gff))

def write_matrix(ofile, matrix):
    """A function to write a dataframe to an output file.
    Args:
        ofile (str) : Path to the output file
        matrix (DataFrame) : a pandas DataFrame
    """
    matrix.to_csv(ofile)



def write_report(ofile, report_data):
    """A function to write the report of the consensus building.
    Args: 
        ofile (str) : Path to the output file
        report_data (dict) : the data to write into the report
    """
    next_section = "-"*24+'\n'
    with open(ofile, 'w') as output_handle:
        output_handle.write(f"#Final number of genes: {report_data['gene_nr']}\n")
        output_handle.write("#Number of genes with support by source:\n")
        for support, number in report_data['support'].items():
            output_handle.write(f"{','.join(support)}\t{number}\n")
        output_handle.write(next_section)
        output_handle.write("#Number of selected genes by source:\n")
        for selected, number in report_data['selected'].items():
            output_handle.write(f"{selected}\t{number}\n")
       

def write_splice(ofile, splice_data, sep='; '):
    """A function to write a splice file for OMA. Used in presence of alternative transcripts for a single gene.
    Args:
        ofile(str) : Path to the file
        splice_data (list) : a list of list, where each of these lists correspond to a gene and contain the protein ID of the product of each of the gene's transcript.
        sep (str) : the separator to use in the splice while, default to semi-clon
    """
    with open(ofile, 'w') as output_handle:
        output_handle.write("\n".join([sep.join(x) for x in splice_data]))

def get_species_tree(codes):
    """A function to obtain a species tree for the exported data from OMA through the OMA Browser API.
    Args:
        codes (list) : a list of OMA's five letter species code. (Ex: HUMAN, MOUSE, DROME)
    Return:
        newick (str) : a newick tree containing the species of interest
    """
    req = urllib3.request('GET', f'https://omabrowser.org/api/taxonomy/?members={",".join(codes)}&type=newick&collapse=true&newick_leaf_label=species_code&newick_internal_label=None&newick_quote_labels=false')
    data = req.data
    newick = json.loads(data)['newick']
    return newick

def edit_parameters_with_species_tree(parameter_file, newick):
    """Modifiy the OMAStandalone parameters file to contain a species tree that contains the exported data. Will need to be manually changed before running next step of OMAnnotation.
    Args:
        parameter_file (str) : Path to the parameter file
        newick (str) : Species tree of exported species, in newick format.
    """
    with open(parameter_file, 'r') as infile:
        all_lines = [line for line in infile.readlines()]
    written_lines = []
    for l in all_lines:
        if l.startswith('SpeciesTree'):
            written_lines.append(f"SpeciesTree := '{newick}'\n")
        else:
            written_lines.append(l)
    with open(parameter_file, 'w') as outfile:
        outfile.write(''.join(written_lines))

def prepare_data(args):
    """Format the input annotations to run OMA Standalone on. First step of OMAnnotation
    Args:
        args : the command line arguments as provided by the user. Expect fasta_folder, gff_annot_folder,  genome_file, splice_folder, OMA_folder and feature_type to have value.
    """
    #Obtain parameters
    input_gff_folder = args.gff_annot_folder
    fasta_folder = args.fasta_folder
    genome_fasta = args.genome_file
    splice_folder = args.splice_folder
    db_folder = os.path.join(args.OMA_folder, 'DB')
    parameter_file = os.path.join(args.OMA_folder,'parameters.drw')
    feature_type_file = args.feature_type
    #Use the correct GFF feature for not standard GFF, extracted from a feature type map provided by the user
    if feature_type_file:
        feature_type_map = extract_feature_map(feature_type_file)
    else:
        feature_type_map = None
    file_list = os.listdir(input_gff_folder)
    spec_codes = []
    #Check the list of species downloaded through OMA
    for pre_oma_data in os.listdir(db_folder):
    #Check for extension (.fa) and a length of 5 - five letter codes
        if pre_oma_data[-3:]=='.fa' and len(pre_oma_data)==8:
            code = pre_oma_data[:-3]
            spec_codes.append(code)
    newick = get_species_tree(spec_codes)
    edit_parameters_with_species_tree(parameter_file, newick)
    #Go through all of the GFF files and make the input files needed for OMA
    for gff in file_list:
        logging.info(f'Starting {gff}')
        prefix = ".".join(gff.split('.')[0:-1])
        logging.info('Reading GFF')
        gff_file = os.path.join(input_gff_folder,gff)
        output_fasta = os.path.join(fasta_folder, prefix+'.fa')
        output_splice= os.path.join(splice_folder, prefix+'.splice')
        fasta_db = os.path.join(db_folder, prefix+'.fa')
        splice_db = os.path.join(db_folder, prefix+'.splice')
        if feature_type_map:
            specific_map = feature_type_map[gff]
        else:
            specific_map = None
        #Read the GFF
        db = extract_gff(gff_file)
        logging.info('Extracting sequences')
        #Extract protein data and splice mapping from the GFF and the genome
        protein_fasta, splice_mapping = extract_peptide_seq(db, genome_fasta, specific_map)
        #Write output files
        logging.info("Writing output file")
        write_fasta(output_fasta, protein_fasta)
        write_fasta(fasta_db, protein_fasta)
        if splice_mapping:
            write_splice(output_splice, splice_mapping)
            write_splice(splice_db, splice_mapping)
        
def extract_feature_map(feature_file):
    """Extract the corresponding feature type in non-standard GFF files (that do not use gene, transcript and CDS feature type)
    The file is user provided and is needed to make sure all GFF are uniformally processed.
    Args:
        feature_file (str) : path to the feature map file
    Returns:
        feature_dict (dict) : an dictionary of standard feature type (str) : feature type of the file (str)
    """
    feature_dict  = dict()
    with open(feature_file,'r') as ffile:
        for line in ffile.readlines():
            line =  line.replace('\n','')
            if line=='':
                continue
            cat = line.split('\t')
            gff_file = cat[0]
            gene_feat = cat[1]
            trans_feat = cat[2]
            cds_feat = cat[3]
            feature_dict[gff_file] = { 'gene' : gene_feat , 'transcript': trans_feat, 'CDS':cds_feat}
    return feature_dict


@timer_wrapper
def extract_peptide_seq(gffdb, fasta, type_map=False):
    """Extract peptide sequence from a genome file using the cds feature from a GFF file by translating the nucleic sequence at the coordinates. This exploit gffutils and Biopython and is based on original code snipets from https://zhiganglu.com/post/py-gff2fasta/
    Args:
        gffdb : a GFF database object from GFF utils
        fasta (str) : path to the genome FASTA file
        type_map (dict) : a dictionary between expected standard GFF features (gene, transcript, CDS) in case the GFF use different nomenclature
    Returns:
        seq_records (list) : a list of BioPython records object for each of the extract protein sequences
        splice_mapping (list) : a list of list, where each of these list correspond to a gene and contain the protein ids of the product of all transcripts of these genes.
    """    
    seq_records = []
    splice_mapping = []
    #Read the genome - use Fasta to read it only once (other .sequence open it anew a very time it is called which slow the software down) 
    fasta = Fasta(fasta)
    #boolean to check if alternative transcripts are predicted in the an notation
    have_splice_variants = False

    if type_map:
        gene_feat = type_map['gene']
    else:
        gene_feat = 'gene'
    #Get gene, transcript and CDS, using type_map if provided
    for g in gffdb.features_of_type(gene_feat, order_by='start'):
        gene = g.id
        transcript_of_gene = []
        if type_map:
            tr_feat = type_map['transcript']
            if type_map['transcript']== type_map['gene']:
                ordered_transcript = [g] 
            else:
                ordered_transcript = list(gffdb.children(g, featuretype=tr_feat, order_by='start'))
        else:
            for transcript_feature in ['mRNA', 'MRNA', 'mrna', 'transcript', 'match']:
                ordered_transcript = list(gffdb.children(g, featuretype=transcript_feature, order_by='start'))
                if len(ordered_transcript)!=0:
                    break

        for t in ordered_transcript:
            transcript_id = t.id
            seq_combined = ''
            #We order the CDS chunk by their starting position to make sure they are in the right order
            if type_map:
                cds_feat = type_map['CDS']
                ordered_child = list(gffdb.children(t, featuretype=cds_feat, order_by='start'))
            else:
                for translated_feature in ['CDS', 'cds']:
                    ordered_child = list(gffdb.children(t, featuretype=translated_feature, order_by='start'))
                    if len(ordered_child)!=0: 
                        break
            #If the strand is negative, CDS are read from the farthest position
            if t.strand=='-':
                ordered_child.reverse()
            start=True
            for i in ordered_child:
                if i.start > i.end:
                    new_end = i.start
                    i.start = i.end
                    i.end = new_end
                #The use_strand option make sure the sequence is read in reverse complement on the reverse strand
                if start:
                    if i.frame=='1':
                        seq_combined+='NN'
                    elif i.frame=='2':
                        seq_combined+='N'
                    start=False
                seq = i.sequence(fasta, use_strand=True) 
                seq_combined += seq
            seq_combined = Seq(seq_combined)
            #Sequence translated using Biopython
            seq_transl = seq_combined.translate()
            #FASTA formatting
            rec = SeqRecord( seq_transl,  id=transcript_id, description=f'gene={gene} transcript={transcript_id}')
            seq_records.append(rec)
            transcript_of_gene.append(transcript_id)
        if len(transcript_of_gene)>1:
            have_splice_variants = True
        splice_mapping.append(transcript_of_gene)
    if not have_splice_variants:
        splice_mapping = None
    return seq_records, splice_mapping


if __name__=='__main__':
    """Main function. Will launch the orchestrator function for either prepare_data or extract_consensus, depending on which one is called (Look above to know more on what they do)
    """
    logging.basicConfig(level=logging.INFO)
    parser = build_arg_parser()
    args = parser.parse_args()
    args.func(args)
