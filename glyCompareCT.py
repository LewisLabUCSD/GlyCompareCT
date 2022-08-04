import argparse
import os
import pandas as pd
import glypy
from glypy.io import glycoct, linear_code, wurcs, iupac
from SPARQLWrapper import SPARQLWrapper
from bs4 import BeautifulSoup
from glycompare import *
import time
import json
from pathlib import PureWindowsPath
import urllib
import subprocess
import requests
import shlex

def main():
    parser = argparse.ArgumentParser(description="GlyCompare command line tool")
    subparsers = parser.add_subparsers(dest="subcommand")
    parser_structure = subparsers.add_parser('structure')
    parser_structure.add_argument("-a", "--abundance", help="Glycan abundance table file", dest="abundance_table", type=str, required=True)
    parser_structure.add_argument("-v", "--var_annot", help="Variable annotation file", dest="variable_annotation", type=str, required=True)
    parser_structure.add_argument("-o", "--output", help="Output directory", dest="output_directory", type=str, required=True)
    parser_structure.add_argument("-s", "--no_linkage", help="Data is pure structural, do not contain linkage information", dest="no_linkage_info", action='store_true')
    parser_structure.add_argument("-c", "--core", help="Number of processors using", dest="num_processors", type=int, default=1)
    parser_structure.add_argument("-p", "--syntax", help="Structural data syntax. Choose from [glycoCT, iupac_extended, linear_code, wurcs, glytoucan_id]", dest="data_syntax", choices=['glycoCT', 'iupac_extended', 'linear_code', 'wurcs', 'glytoucan_id'], required=True)
    parser_structure.add_argument("-n", "--norm", help="Glycan abundance normalization. Choose from [none, min-max, prob-quot]", dest="glycan_abundance_normalization", choices=['none', 'min-max', 'prob-quot'], default='none')
    parser_structure.add_argument("-m", "--multiplier", help="Substructure abundance multiplier. Choose from [binary, integer]", dest="multiplier", choices=['binary', 'integer'], default='integer')
    parser_structure.add_argument("-b", "--no_sub_norm", help="Do not normalize substructure abundance. Leave as the absolute value", dest="no_substructure_normlization", action='store_true')
    # custom root only used when root is set to custom
    parser_structure.add_argument("-r", "--root", help="Set the glycan root. Choose from [epitope, N, O, lactose, custom]", dest="root", choices=['epitope', 'N', 'O', 'lactose', 'custom'], default='epitope')
    parser_structure.add_argument("-u", "--custom_root", help="Custom root", dest="custom_root", type=str, default='')
    parser_structure.add_argument("-d", "--heatmap", help="Draw motif abundance heatmap", dest="heatmap", action='store_true')
    parser_structure.add_argument("-i", "--ignore", help="Ignore non-recognized glycan structures and proceed the rest", dest="ignore", action='store_true')
    parser_structure.add_argument("-e", "--share", help="Either run locally or register the output motif structures to Glytoucan. Choose from [private, register]", dest="share", choices=["private", "register"], default="private")
    parser_structure.add_argument("-C", "--Contributor_ID", help="Glytoucan Contributor ID", dest="CID", type=str)
    parser_structure.add_argument("-A", "--API_key", help="Glytoucan API key", dest="API", type=str)
    
    
    parser_composition = subparsers.add_parser('composition')
    parser_composition.add_argument("-a", "--abundance", help="Glycan abundance table file", dest="abundance_table", type=str, required=True)
    parser_composition.add_argument("-v", "--var_annot", help="Variable annotation file", dest="variable_annotation", type=str, required=True)
    parser_composition.add_argument("-o", "--output", help="Output directory", dest="output_directory", type=str, required=True)
    parser_composition.add_argument("-n", "--norm", help="Glycan abundance normalization. Choose from [none, min-max, prob-quot]", dest="glycan_abundance_normalization", choices=['none', 'min-max', 'prob-quot'], default='none')
    parser_composition.add_argument("-i", "--ignore", help="Ignore non-recognized glycan compositions and proceed the rest", dest="ignore", action='store_true')

    parser_annotate = subparsers.add_parser('annotate')
    parser_annotate.add_argument("-n", "--annotation_table", help="Motif annotation table file", dest="annotation_table", type=str, required=True)
    
    args = parser.parse_args()
    if not args.subcommand:
        parser.parse_args(["--help"])
    elif args.subcommand == 'structure':
        input_validation(args)
        parser.set_defaults(func=structure)
        args = parser.parse_args()
        args.func(args)
    elif args.subcommand == 'composition':
        input_validation(args)
        parser.set_defaults(func=composition)
        args = parser.parse_args()
        args.func(args)
    elif args.subcommand == 'annotate':
        input_validation(args)
        parser.set_defaults(func=annotation)
        args = parser.parse_args()
        args.func(args)

# Validate input:
# abundance table: 0. file path 1. columns uniqueness (glycans) 2. rows uniqueness (samples) 2. non-negativity
# variable annotation: 0. file path 1. columns (Name + Glycan Structure) 2. structure data format (whether it's the input syntax) or composition data format
def input_validation(args):

    if args.subcommand == "annotate":
        annot = pd.read_csv(args.annotation_table)
        # Validate correct format of annotation file
        assert "WURCS" in annot.columns or "glycoCT" in annot.columns, "Annotation table needs to have at least one of WURCS or glycoCT column."
        return
    
    # Validate custom root is set if root is set to be custom.
    if args.subcommand == "structure" and args.root == 'custom' and not args.custom_root:
        raise Exception("Please specify custom root -u")

    if args.subcommand == "structure" and args.share == 'register' and (not args.CID or not args.API):
        raise Exception("You need to set Glytoucan contributor ID and API key to use register mode. You can get them by signing up on https://glytoucan.org/")
        
    print("Validating input files...")
    ### Validate file paths

    # Validate glycan abundance table file path
    assert os.path.isfile(args.abundance_table), "Invalid path to glycan abundance table csv file, check the path: " + args.abundance_table

    # Validate variable annotation file path
    assert os.path.isfile(args.variable_annotation), "Invalid path to variable annotation csv file, check the path: " + args.variable_annotation

    # Validate output file path
    assert os.path.isdir(os.sep.join(args.output_directory.split(os.sep)[:-1])), "Invalid output path, check the path: " + os.sep.join(args.output_directory.split(os.sep)[:-1])

    # Validate file names
    assert ".csv" in args.abundance_table and ".csv" in args.variable_annotation, "Invalid file type. Only csv files can be used as abundance table and variable annotation"

#     assert "_abundance_table.csv" in args.abundance_table and "_variable_annotation.csv" in args.variable_annotation and args.abundance_table.split("/")[-1].split("_abundance_table.csv")[0] == args.variable_annotation.split("/")[-1].split("_variable_annotation.csv")[0], "Invalid naming of input files. The files must be .csv files and their prefix should match. The suffix should be _abundance_table.csv and _variable_annotation.csv respectively"

    ### Validate glycan abundance table

    # Validate column glycan names uniqueness
    f = open(args.abundance_table, "r")
    col = f.readlines()[0].split("\n")[0]
    assert len(set(col.split(","))) == len(col.split(",")), "Duplicate column glycan names found: " + str(set([i for i in col.split(",") if i and col.split(",").count(i) != 1]))

    # Validate row sample names uniqueness
    glycan_abd = pd.read_csv(args.abundance_table, index_col = 0)
    assert len(set(list(glycan_abd.index))) == len(list(glycan_abd.index)), "Duplicate row sample names found: " + str(set([i for i in glycan_abd.index if list(glycan_abd.index).count(i) != 1]))

    # Validate non-negativity of abundance
    assert (glycan_abd < 0).sum().sum() == 0, "The input glycan abundance table can only have non-negative values. Negative values detected. Please correct."

    ### Validate variable annotation

    var_annot = pd.read_csv(args.variable_annotation)
    # Validate structural data
    if args.subcommand == 'structure':
        # Validate column names
        assert "Name" in var_annot.columns, "'Name' not found in variable annotation column names"
        assert "Glycan Structure" in var_annot.columns, "'Glycan Structure' not found in variable annotation column names"

        # Validate the Name column of the variable annotation file is the same as the column names of glycan abundance table
        assert set(glycan_abd.columns) == set(var_annot["Name"]), "The 'Name' column values of the variable annotation file are not the same as the column names of the glycan abundance table. Please check and correct it."

        # Validate the data format
        bad_glycan_names = []
        for i in range(len(var_annot["Glycan Structure"])):
            if not glycan_syntax_validation(var_annot["Glycan Structure"].tolist()[i], args.data_syntax):
                bad_glycan_names.append(var_annot["Name"].tolist()[i])
        if bad_glycan_names:
            output_path = os.path.abspath(os.sep.join(args.output_directory.split(os.sep)[:-1]))
            project_name = args.output_directory.split(os.sep)[-1] if args.output_directory.split(os.sep)[-1] else "glyCompareCT"
            output_path = os.path.abspath(os.path.join(output_path, project_name))
            if not os.path.exists(output_path):
                os.makedirs(output_path)
            f = open(output_path + os.sep + "bad_glycans.txt", "w")
            f.write(",".join(bad_glycan_names))
            f.close()
            if not args.ignore:
                assert not bad_glycan_names, ",".join(bad_glycan_names) + "\n" + str(len(bad_glycan_names)) + " of " + str(var_annot.shape[0]) + " glycans have invalid gstructure syntax. The glycan structures of the above names failed to be recognized as " + str(args.data_syntax) + ". Non-recognized glycans are saved to bad_glycans.txt. Consider using -i to ignore non-recognized glycans and proceed."
            else:
                print("\x1b[33;20m Warning: Ignored "  + str(len(bad_glycan_names)) + " of " + str(var_annot.shape[0]) +  " glycans that have invalid structure syntax. \033[0m")
                if not os.path.exists(output_path):
                    os.makedirs(output_path)
                glycan_abd.drop(bad_glycan_names, axis = 1).to_csv(output_path + os.sep + "temp_abundance_table.csv", index=True)
                var_annot.drop(var_annot[var_annot["Name"].isin(bad_glycan_names)].index, axis = 0).to_csv(output_path + os.sep + "temp_variable_annotation.csv", index=False)

        # Validate optional custom root file
        if args.custom_root:
            assert os.path.isfile(args.custom_root), "Invalid path to custom root file, check the path: " + args.custom_root
            f = open(args.custom_root, "r")
            root = "".join(f.readlines())
            try:
                glycoct.loads(root)
            except:
                raise Exception("Invalid custom root: " + root)

    # Validate compositional data
    elif args.subcommand == 'composition':
        # Validate column names
        assert "Name" in var_annot.columns, "'Name' not found in variable annotation column names"
        assert "Composition" in var_annot.columns, "'Composition' not found in variable annotation column names"

        # Validate the Name column of the variable annotation file is the same as the column names of glycan abundance table
        assert set(glycan_abd.columns) == set(var_annot["Name"]), "The 'Name' column values of the variable annotation file are not the same as the column names of the glycan abundance table. Please check and correct it."

        # Validate the data format
#         iupac_syms = ['GlcNAc', 'GalNAc', 'HexNAc', 'Neu5Ac', 'GalA', 'GlcA', 'Glc', 'Gal', 'Man', 'Neu', 'NAc', 'KDN', 'Kdo', 'Ido', 'Rha', 'Fuc', 'Xyl', 'Rib', 'Ara', 'All', 'Api', 'Fru', 'Hex', '4eLeg', '6dAlt', '6dAltNAc', '6dGul', '6dTal', '6dTalNAc', 'Abe', 'Aci', 'AllA', 'AllN', 'AllNAc', 'Alt', 'AltA', 'AltN', 'AltNAc', 'Bac', 'Col', 'DDmanHep', 'Dha', 'Dig', 'FucNAc', 'GalN', 'GlcN', 'Gul', 'GulA', 'GulN', 'GulNAc', 'IdoA', 'IdoN', 'IdoNAc', 'Kdn', 'Leg', 'LDmanHep', 'Lyx', 'ManA', 'ManN', 'ManNAc', 'Mur', 'MurNAc', 'MurNGc', 'Neu5Gc', 'Oli', 'Par', 'Pse', 'Psi', 'Qui', 'QuiNAc', 'RhaNAc', 'Sia', 'Sor', 'Tag', 'Tal', 'TalA', 'TalN', 'TalNAc', 'Tyv', 'Phospho', 'NeuAc']
        bad_glycan_names = []
        for i in range(len(var_annot["Composition"])):
            if not glycan_syntax_validation(var_annot["Composition"].tolist()[i], "composition"):
                bad_glycan_names.append(var_annot["Name"].tolist()[i])
        if bad_glycan_names:
            output_path = os.path.abspath(os.sep.join(args.output_directory.split(os.sep)[:-1]))
            project_name = args.output_directory.split(os.sep)[-1] if args.output_directory.split(os.sep)[-1] else "glyCompareCT"
            output_path = os.path.abspath(os.path.join(output_path, project_name))
            if not os.path.exists(output_path):
                os.makedirs(output_path)
            f = open(output_path + os.sep + "bad_glycans.txt", "w")
            f.write(",".join(bad_glycan_names))
            f.close()
            if not args.ignore:
                assert not bad_glycan_names, str(len(bad_glycan_names)) + " of " + str(var_annot.shape[0]) + " glycans have invalid glycan composition format. The following glycan compositions failed to be recognized as composition data: \n" + "\,".join(bad_glycan_names) + "\nPlease check monosaccharide names and parenthesis closure. Non-recognized glycans are saved to bad_glycans.txt. Consider using -i to ignore non-recognized glycans and proceed."
            else:
                print("\x1b[33;20m Warning: Ignored "  + str(len(bad_glycan_names)) + " of " + str(var_annot.shape[0]) +  " glycans that have invalid composition syntax. \033[0m")
                if not os.path.exists(output_path):
                    os.makedirs(output_path)
                glycan_abd.drop(bad_glycan_names, axis = 1).to_csv(output_path + os.sep + "temp_abundance_table.csv", index=True)
                var_annot.drop(var_annot[var_annot["Name"].isin(bad_glycan_names)].index, axis = 0).to_csv(output_path + os.sep + "temp_variable_annotation.csv", index=False)
            

# Fetch glytoucan ID in annotation table. 
def annotation(args):
    print("Start annotation mode")
    annot_path = os.path.abspath(args.annotation_table)
    annot = pd.read_csv(annot_path, index_col = 0)
    new_accs = []
    accs = []
    print(str(annot.shape[0]) + " glycans in total")
    if "Glytoucan ID" in annot.columns:
        accs = annot["Glytoucan ID"].tolist()
    print("Annotation started")
    for i in range(annot.shape[0]):
        if i % 5 == 0 and i != 0:
            print(str(len(new_accs) - new_accs.count("na")) + "/" + str(i) + " glytoucan ID annotated")
        if accs and accs[i] != "na":
            new_accs.append(accs[i])
        elif "WURCS" in annot.columns:
            
            new_accs.append(fetch_glytoucan(annot["WURCS"].tolist()[i], "WURCS"))
        elif "glycoCT" in annot.columns:
            new_accs.append(fetch_glytoucan(annot["glycoCT"].tolist()[i], "glycoCT"))
    print(str(len(new_accs) - new_accs.count("na")) + "/" + str(annot.shape[0]) + " glytoucan ID annotated in total")
    annot["Glytoucan ID"] = new_accs
    annot.to_csv(annot_path)
    
    

# Running glyCompare on structural data
def structure(args):
    print("Start structure mode")
    print("Initializing GlyCompare...")
    abd_path = os.path.abspath(args.abundance_table)
    var_path = os.path.abspath(args.variable_annotation)
    output_path = os.path.abspath(os.sep.join(args.output_directory.split(os.sep)[:-1]))
    project_name = args.output_directory.split(os.sep)[-1] if args.output_directory.split(os.sep)[-1] else "glyCompareCT"
    output_path = os.path.abspath(os.path.join(output_path, project_name))
    if os.path.isfile(output_path + os.sep + "temp_variable_annotation.csv") and os.path.isfile(output_path + os.sep + "temp_abundance_table.csv") and args.ignore:
        abd_path = output_path + os.sep + "temp_abundance_table.csv"
        var_path = output_path + os.sep + "temp_variable_annotation.csv"
    working_addr = output_path
    reference_addr = os.path.join(os.path.dirname(os.path.abspath(__file__)), "reference")
    keywords_dict = pipeline_functions.load_para_keywords(project_name, working_addr, reference_addr = reference_addr)
    # Create temperary source data. Delete after GlyCompare is done.
    pipeline_functions.check_init_dir(keywords_dict)
    if os.sep == "/":
        os.popen("cp " + "\ ".join(abd_path.split(" ")) + " " + "\ ".join(output_path.split(" ")) + "/source_data/" + project_name + "_abundance_table.csv")
        os.popen("cp " + "\ ".join(var_path.split(" ")) + " " + "\ ".join(output_path.split(" ")) + "/source_data/" + project_name + "_variable_annotation.csv")
    else:
        os.popen("copy \"" + abd_path + "\" \"" + output_path + "\\source_data\\" + project_name + "_abundance_table.csv" + "\"")
        os.popen("copy \"" + var_path + "\" \"" + output_path + "\\source_data\\" + project_name + "_variable_annotation.csv" + "\"")

    while True:
        time.sleep(5)
        if (os.path.isfile(output_path + "/source_data/" + project_name + "_abundance_table.csv") and os.path.isfile(output_path + "/source_data/" + project_name + "_variable_annotation.csv")) or (os.path.isfile(output_path + "\\source_data\\" + project_name + "_abundance_table.csv") and os.path.isfile(output_path + "\\source_data\\" + project_name + "_variable_annotation.csv")):
            break


    print("Generating glycoCT local files and glycan dictionary...")
    syntax_convert = {"glycoCT": "glycoCT", "iupac_extended": "IUPAC_extended", "linear_code": "linear code", "wurcs": "WURCS", "glytoucan_id": "glytoucan ID"}
    glycan_type = syntax_convert[args.data_syntax]
    pipeline_functions.generate_glycoct_files(keywords_dict, glycan_type)
    if args.no_linkage_info:
        linkage_specific = False
        merged_list = [keywords_dict['structure_only_glycoct_reference'], keywords_dict['structure_only_wurcs_reference']]
        reference_dict_addr = keywords_dict['structure_only_reference']
    else:
        linkage_specific = True
        merged_list = [keywords_dict['linkage_specific_glycoct_reference'], keywords_dict['linkage_specific_wurcs_reference']]
        reference_dict_addr = keywords_dict['linkage_specific_reference']
    reference_dict = json.load(open(reference_dict_addr, "r"))
    reverse_dict = {v: k for k, v in reference_dict.items()}
    glycan_dict = pipeline_functions.load_glycans_pip(keywords_dict = keywords_dict, data_type='local_glycoct', reference_dict = reference_dict, reverse_dict = reverse_dict, linkage_specific = linkage_specific, reference_dict_addr = reference_dict_addr)

    print("Creating glycan_substructure_occurance_dict...")
    keywords_dict = pipeline_functions.load_para_keywords(project_name, working_addr, reference_addr = reference_addr)
    if args.no_linkage_info:
        reference_dict_addr = keywords_dict['structure_only_reference']
    else:
        reference_dict_addr = keywords_dict['linkage_specific_reference']
    reference_dict = json.load(open(reference_dict_addr, "r"))
    reverse_dict = {v: k for k, v in reference_dict.items()}
    matched_df = pipeline_functions.extract_and_merge_substrutures_pip(keywords_dict, num_processors=args.num_processors, reference_dict = reference_dict, reverse_dict = reverse_dict, linkage_specific=linkage_specific, forced=True, merged_list = merged_list, reference_dict_addr = reference_dict_addr)
    matched_df = ""
    keywords_dict = pipeline_functions.load_para_keywords(project_name, working_addr, reference_addr = reference_addr)
    if args.no_linkage_info:
        reference_dict_addr = keywords_dict['structure_only_reference']
    else:
        reference_dict_addr = keywords_dict['linkage_specific_reference']

    reference_dict = json.load(open(reference_dict_addr, "r"))
    reverse_dict = {v: k for k, v in reference_dict.items()}
    chra_to_id = {}
    for i in glycan_dict.keys():
        chra_to_id[i]=i
    json_utility.store_json(keywords_dict['name_to_id_addr'], chra_to_id)
    glycan_dict = ""

    print("Creating glycoprofile_list...")
    glycan_abd_table = pd.read_csv(os.path.join(keywords_dict['source_dir'], project_name + '_abundance_table.csv'), index_col=0)
    glycan_abd_table = glycan_abd_table[(glycan_abd_table.T != 0).any()]
    norm = args.glycan_abundance_normalization
    if norm == "min-max":
        glycan_abd_table = pipeline_functions.normalization(glycan_abd_table, style = "std")
    elif norm == "prob-quot":
        glycan_abd_table = pipeline_functions.normalization(glycan_abd_table, style = "pq")
    elif norm == "none":
        glycan_abd_table = glycan_abd_table
    glycan_abd_table = glycan_abd_table.transpose()
    pd.set_option('use_inf_as_na', True)
    glycan_abd_table = glycan_abd_table.dropna(axis=0, how="all")
    glycan_abd_table = glycan_abd_table.dropna(axis=1, how="all")
    if args.no_substructure_normlization:
        absolute = True
    else:
        absolute = False
    if args.multiplier == "binary":
        get_existance = True
    elif args.multiplier == "integer":
        get_existance = False
    _, glycoprofile_list = pipeline_functions.glycoprofile_pip(keywords_dict, glycan_abd_table, unique_glycan_identifier_to_structure_id=True, already_glytoucan_id=False, external_profile_naming=True, forced=True, absolute = absolute, get_existance = get_existance)
    glycan_abd_table = ""
    glycoprofile_list = ""

    print("Creating motif abundance table...")
    core_input = args.root
    if core_input == "N":
        core = glycoct.loads('RES\n1b:b-dglc-HEX-1:5\n2s:n-acetyl\n3b:b-dglc-HEX-1:5\n4s:n-acetyl\n5b:b-dman-HEX-1:5\nLIN\n1:1d(2+1)2n\n2:1o(4+1)3d\n3:3d(2+1)4n\n4:3o(4+1)5d\n')
        only_substructures_start_from_root = True
    elif core_input == "O":
        core = glycoct.loads("""
        RES
        1b:b-dgal-HEX-x:x
        2s:n-acetyl
        LIN
        1:1d(2+1)2n""")
        only_substructures_start_from_root = True
    elif core_input == "lactose":
        core = glycoct.loads("""
        RES
        1b:b-dglc-HEX-1:5
        2b:b-dgal-HEX-1:5
        LIN
        1:1o(4+1)2d""")
        only_substructures_start_from_root = True
    elif core_input == "custom":
        f = open(args.custom_root, "r")
        temp = "".join(f.readlines())
        core = glycoct.loads(temp)
        only_substructures_start_from_root = True
    elif core_input == "epitope":
        core = ""
        only_substructures_start_from_root = False
    motif_abd_table, motif_lab, merged_weights_dict = pipeline_functions.select_motifs_pip(keywords_dict, linkage_specific=linkage_specific, only_substructures_start_from_root=only_substructures_start_from_root, reverse_dict = reverse_dict, num_processors=args.num_processors, core=core, drop_parellel=True, drop_diff_abund=False, select_col= [])

    reference_vector = json.load(open(reference_dict_addr, "r"))
    motif_glycoct = json.load(open(keywords_dict['motif_glycoct_dict_addr'], "r"))
    motif_names = {}
    for key in motif_glycoct.keys():
        ref_name = motif_glycoct[key]
        motif_names[key] = ref_name
    index_col = list(motif_abd_table.index)
    motif_abd_table.index = [motif_names[str(i)] for i in list(motif_abd_table.index)]
    ref_col = list(motif_abd_table.index)
    
    glycoct_col = [reverse_dict[motif_glycoct[str(i)]] for i in index_col]
    wurcs_col = [wurcs.dumps(glycoct.loads(i)) for i in glycoct_col]
    motif_abd_table_addr = keywords_dict['motif_abd_table_addr']
    motif_annotation_addr = motif_abd_table_addr.split(project_name + "_")[0] + project_name + "_motif_annotation.csv"
    motif_abd_table.to_csv(motif_abd_table_addr)
    if args.share == "register":
        print("Fetching glytoucan ID...")
        glytoucan_col = []
        print(str(len(wurcs_col)) + " glycans in total")
        for i in range(len(wurcs_col)):
            if i % 5 == 0 and i != 0:
                print(str(len(glytoucan_col) - glytoucan_col.count("na")) + "/" + str(i) + " glytoucan ID annotated")
            glytoucan_col.append(fetch_glytoucan(wurcs_col[i], "WURCS"))
        print(str(len(glytoucan_col) - glytoucan_col.count("na")) + "/" + str(len(wurcs_col)) + " glytoucan ID annotated in total")
        
        motif_annotation = pd.DataFrame(zip(index_col, glytoucan_col, wurcs_col, glycoct_col), columns = ["Substructure Index", "Glytoucan ID", "WURCS", "glycoCT"])
        motif_annotation.index = ref_col
    else:
        motif_annotation = pd.DataFrame(zip(index_col, wurcs_col, glycoct_col), columns = ["Substructure Index", "WURCS", "glycoCT"])
        motif_annotation.index = ref_col
    
    
    motif_annotation.to_csv(motif_annotation_addr)


    if args.heatmap:
        print("Drawing heatmap...")
        motif_abd_table = motif_abd_table.dropna(axis = 1)
        index = list(motif_abd_table.index)
        unique_rows = motif_abd_table.stack().groupby(level=0).apply(lambda x: x.unique().tolist())
        for i in range(len(index)):
            if len(unique_rows[index[i]]) == 1:
                print("Row " +  str(index[i]) + " is dropped from motif_abd_table becuase it contains same values")
                motif_abd_table = motif_abd_table.drop([index[i]])
        glycoprofile_cluster_dict, glyco_motif_cluster_dict = pipeline_functions.clustering_analysis_pip(keywords_dict=keywords_dict, motif_abd_table=motif_abd_table, select_profile_name=[])
        os.rename(keywords_dict['plot_output_dir'] + "pseudo_profile_clustering.svg", keywords_dict['plot_output_dir'] + project_name + "_pseudo_profile_clustering.svg")
        os.rename(keywords_dict['plot_output_dir'] + "profile_clustering.svg", keywords_dict['plot_output_dir'] + project_name + "_profile_clustering.svg")
        os.rename(keywords_dict['plot_output_dir'] + "motif_cluster.svg", keywords_dict['plot_output_dir'] + project_name + "_motif_cluster.svg")

    print("Removing intermediate data...")
    pre_dt = output_path + os.sep + "glycoct" + os.sep
    for file in os.listdir(pre_dt):
        if os.path.isfile(pre_dt + file):
            os.remove(pre_dt + file)
    os.rmdir(pre_dt)
    pre_sr = output_path + os.sep + "source_data" + os.sep
    for file in os.listdir(pre_sr):
        if os.path.isfile(pre_sr + file):
            os.remove(pre_sr + file)
    os.rmdir(pre_sr)

    if os.path.isdir(output_path + os.sep + project_name + "_output_data"):
        pre_ot = output_path + os.sep + project_name + "_output_data" + os.sep
        for file in os.listdir(pre_ot):
            if os.path.isfile(pre_ot + file):
                os.remove(pre_ot + file)
        os.rmdir(pre_ot)
    os.rename(output_path + os.sep + "output_data", output_path + os.sep + project_name + "_output_data")
    
    if args.share == "register":
#         pp = os.path.abspath(output_path + os.sep + project_name + "_output_data/" + project_name + "_register_submission.txt")
        pp = "." + os.sep + project_name + "_register_submission.txt"
#         if os.sep == "/":
#             pp = pp.replace(' ', '\ ')
#         else:
#             pp = "\"" + pp + "\""
        f = open(pp, "w")
        rc = 0
        for i in range(len(glytoucan_col)):
            if glytoucan_col[i] == "na":
                f.write(wurcs_col[i] + "\n")
                rc += 1
        f.close()
        
        if rc:
            print("Registering glycans...")
#             command = ["bash", "." + os.sep + "curlRegisterFile.sh"]
            command = "bash" + " curlRegisterFile.sh " + args.CID + " " + args.API
#             print(command.split() + [pp])
            process = subprocess.Popen(command.split() + [pp], stdout=subprocess.PIPE, text=True)
            output, error = process.communicate()
#             command_result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
#             process.stdin.write(args.CID + " " + args.API + " " + pp)
#             process.stdin.close()
#             if process.stderr.read():
#                 print(process.stderr.read())
#             print(process.stdout.read())
            if error:
                print(error)
            print(output)
        os.remove(pp)

    if not any(os.scandir(output_path + os.sep + "output_plot")):
        os.rmdir(output_path + os.sep + "output_plot")
    else:
        if os.path.isdir(output_path + os.sep + project_name + "_output_plot"):
            pre_pl = output_path + os.sep + project_name + "_output_plot" + os.sep
            for file in os.listdir(pre_pl):
                if os.path.isfile(pre_pl + file):
                    os.remove(pre_pl + file)
            os.rmdir(pre_pl)
        os.rename(output_path + os.sep + "output_plot", output_path + os.sep + project_name + "_output_plot")
    
    if os.path.isfile(output_path + os.sep + "temp_abundance_table.csv"):
        os.remove(output_path + os.sep + "temp_abundance_table.csv")
    if os.path.isfile(output_path + os.sep + "temp_variable_annotation.csv"):
        os.remove(output_path + os.sep + "temp_variable_annotation.csv")


# Running glyCompare on compositional data
def composition(args):
    print("Start composition mode...")
    print("Initializing GlyCompare...")
    abd_path = os.path.abspath(args.abundance_table)
    var_path = os.path.abspath(args.variable_annotation)
    output_path = os.path.abspath(os.sep.join(args.output_directory.split(os.sep)[:-1]))
    project_name = args.output_directory.split(os.sep)[-1] if args.output_directory.split(os.sep)[-1] else "glyCompareCT"
    output_path = os.path.abspath(os.path.join(output_path, project_name))
    if os.path.isfile(output_path + os.sep + "temp_variable_annotation.csv") and os.path.isfile(output_path + os.sep + "temp_abundance_table.csv") and args.ignore:
        abd_path = output_path + os.sep + "temp_abundance_table.csv"
        var_path = output_path + os.sep + "temp_variable_annotation.csv"
    working_addr = output_path

    reference_addr = os.path.join(os.path.dirname(os.path.abspath(__file__)), "reference")
    keywords_dict = pipeline_functions.load_para_keywords(project_name, working_addr, reference_addr)
    pipeline_functions.check_init_dir(keywords_dict)
    # Create temperary source data. Delete after GlyCompare is done.
    if os.sep == "/":
        os.popen("cp " + "\ ".join(abd_path.split(" ")) + " " + "\ ".join(output_path.split(" ")) + "/source_data/" + project_name + "_abundance_table.csv")
        os.popen("cp " + "\ ".join(var_path.split(" ")) + " " + "\ ".join(output_path.split(" ")) + "/source_data/" + project_name + "_variable_annotation.csv")
    else:
        os.popen("copy \"" + abd_path + "\" \"" + output_path + "\\source_data\\" + project_name + "_abundance_table.csv" + "\"")
        os.popen("copy \"" + var_path + "\" \"" + output_path + "\\source_data\\" + project_name + "_variable_annotation.csv" + "\"")

    while True:
        time.sleep(5)
        if (os.path.isfile(output_path + "/source_data/" + project_name + "_abundance_table.csv") and os.path.isfile(output_path + "/source_data/" + project_name + "_variable_annotation.csv")) or (os.path.isfile(output_path + "\\source_data\\" + project_name + "_abundance_table.csv") and os.path.isfile(output_path + "\\source_data\\" + project_name + "_variable_annotation.csv")):
            break

    print("Creating motif abundance table...")
    protein_sites = "all"
    norm = args.glycan_abundance_normalization
    if norm == "min-max":
        norm_parsed = "min-max"
    elif norm == "prob-quot":
        norm_parsed = "pq"
    elif norm == "none":
        norm_parsed = "no"
    motif_abd, directed_edge_list = pipeline_functions.compositional_data(keywords_dict, protein_sites = protein_sites, reference_vector = None, forced = True, norm = norm_parsed)
    pre_dt = output_path + os.sep + "glycoct" + os.sep
    for file in os.listdir(pre_dt):
        if os.path.isfile(pre_dt + file):
            os.remove(pre_dt + file)
    os.rmdir(pre_dt)
    pre_sr = output_path + os.sep + "source_data" + os.sep
    for file in os.listdir(pre_sr):
        if os.path.isfile(pre_sr + file):
            os.remove(pre_sr + file)
    os.rmdir(pre_sr)

    if os.path.isdir(output_path + os.sep + project_name + "_output_data"):
        pre_ot = output_path + os.sep + project_name + "_output_data" + os.sep
        for file in os.listdir(pre_ot):
            if os.path.isfile(pre_ot + file):
                os.remove(pre_ot + file)
        os.rmdir(pre_ot)
    os.rename(output_path + os.sep + "output_data", output_path + os.sep + project_name + "_output_data")

    if not any(os.scandir(output_path + os.sep + "output_plot")):
        os.rmdir(output_path + os.sep + "output_plot")
    else:
        if os.path.isdir(output_path + os.sep + project_name + "_output_plot"):
            pre_pl = output_path + os.sep + project_name + "_output_plot" + os.sep
            for file in os.listdir(pre_pl):
                if os.path.isfile(pre_pl + file):
                    os.remove(pre_pl + file)
            os.rmdir(pre_pl)
        os.rename(output_path + os.sep + "output_plot", output_path + os.sep + project_name + "_output_plot")


# Validate glycan data syntax.
# syntax: glycoCT, iupac_extended, linear_code, wurcs, glytoucan_id
def glycan_syntax_validation(data, syntax):
    if syntax == "glycoCT":
        try:
            temp = glycoct.loads(data)
            glycoct.dumps(temp)
            assert isinstance(temp, glypy.Glycan)
        except:
            return False
    elif syntax == "iupac_extended":
        try:
            temp = iupac.loads(data)
            iupac.dumps(temp)
            assert isinstance(temp, glypy.Glycan)
        except:
            return False
    elif syntax == "linear_code":
        try:
            temp = linear_code.loads(data)
            linear_code.dumps(temp)
            assert isinstance(temp, glypy.Glycan)
        except:
            return False
    elif syntax == "wurcs":
        try:
            temp = wurcs.loads(data)
            wurcs.dumps(temp)
            assert isinstance(temp, glypy.Glycan)
        except:
            return False
    # Seems like glytoucan SPARQL query is not updated as glytoucan web app. So latest data are likely not querable.
    elif syntax == "glytoucan_id":
        try:
            gct = get_glycoct_from_glytoucan(data)
            glycoct.dumps(glycoct.loads(gct))
            assert isinstance(glycoct.loads(gct), glypy.Glycan)
        except:
            return False
    elif syntax == "composition":
#         iupac_syms = ['GlcNAc', 'GalNAc', 'HexNAc', 'Neu5Ac', 'GalA', 'GlcA', 'Glc', 'Gal', 'Man', 'Neu', 'NAc', 'KDN', 'Kdo', 'Ido', 'Rha', 'Fuc', 'Xyl', 'Rib', 'Ara', 'All', 'Api', 'Fru', 'Hex', '4eLeg', '6dAlt', '6dAltNAc', '6dGul', '6dTal', '6dTalNAc', 'Abe', 'Aci', 'AllA', 'AllN', 'AllNAc', 'Alt', 'AltA', 'AltN', 'AltNAc', 'Bac', 'Col', 'DDmanHep', 'Dha', 'Dig', 'FucNAc', 'GalN', 'GlcN', 'Gul', 'GulA', 'GulN', 'GulNAc', 'IdoA', 'IdoN', 'IdoNAc', 'Kdn', 'Leg', 'LDmanHep', 'Lyx', 'ManA', 'ManN', 'ManNAc', 'Mur', 'MurNAc', 'MurNGc', 'Neu5Gc', 'Oli', 'Par', 'Pse', 'Psi', 'Qui', 'QuiNAc', 'RhaNAc', 'Sia', 'Sor', 'Tag', 'Tal', 'TalA', 'TalN', 'TalNAc', 'Tyv', 'Phospho', 'NeuAc']
        decomp = [i.split("(") for i in data.split(")")]
        if not (decomp[-1] == [''] and sum([type(int(i[1])) is not int for i in decomp[:-1]]) == 0):
            return False
    else:
        raise Exception("Invalid syntax " + syntax + ". Choose among glycoCT, iupac_extended, linear_code, wurcs, and glytoucan_id")
    return True


def fetch_glytoucan(gly, syn):
    if syn == "WURCS":
        tt = urllib.parse.quote_plus(gly)
    elif syn == "glycoCT":
        ww = wurcs.dumps(glycoct.loads(gly))
        tt = urllib.parse.quote_plus(ww)
    r = requests.get("https://api.glycosmos.org/sparqlist/wurcs2gtcids?wurcs=%s" % tt)
    soup = BeautifulSoup(r.text, 'html.parser')
    if "id" not in soup.text or "wurcs" not in soup.text:
        acc = "na"
    else:
        acc = soup.text.split(": \"")[1].split("\",")[0]
    return acc
    
def get_wurcs_from_glytoucan(ID):
    sparql = SPARQLWrapper("https://ts.glytoucan.org/sparql")
    sparql.setQuery("""
    PREFIX glycan: <http://purl.jp/bio/12/glyco/glycan#>
    PREFIX glytoucan: <http://www.glytoucan.org/glyco/owl/glytoucan#>
    SELECT DISTINCT ?WURCS_label
    WHERE {
      # Accession Number
      ?saccharide glytoucan:has_primary_id ?accNum .
      FILTER (?accNum = '%s')
      # WURCS
      OPTIONAL{
      ?saccharide glycan:has_glycosequence ?wcsSeq .
      ?wcsSeq glycan:has_sequence ?wcsLabel .
      BIND(STR(?wcsLabel) AS ?WURCS_label)
      ?wcsSeq glycan:in_carbohydrate_format glycan:carbohydrate_format_wurcs .
      }
    }
    """ % ID)

    try:
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        structure = results["results"]["bindings"][0]['WURCS_label']['value']
    except:
        print("No WURCS structure found for ", ID)
        structure = ''
    return structure

def get_glycoct_from_glytoucan(ID):
    # Returns glycan GlycoCT structure from glytoucan ID
    # Perform the query
    AccNum = '"' + ID + '"'
    sparql = SPARQLWrapper("http://ts.glytoucan.org/sparql")
    query = '''
PREFIX glycan: <http://purl.jp/bio/12/glyco/glycan#>
PREFIX glytoucan:  <http://www.glytoucan.org/glyco/owl/glytoucan#>
SELECT DISTINCT ?Sequence
FROM <http://rdf.glytoucan.org/core>
FROM <http://rdf.glytoucan.org/sequence/glycoct>
WHERE {
    VALUES ?PrimaryId {%s}
    ?Saccharide glytoucan:has_primary_id ?PrimaryId .
    ?Saccharide glycan:has_glycosequence ?GlycoSequence .
    ?GlycoSequence glycan:has_sequence ?Sequence .
    ?GlycoSequence glycan:in_carbohydrate_format glycan:carbohydrate_format_glycoct.
    }
    ''' % AccNum
    sparql.setQuery(query)
    results = sparql.query().convert()

    # Parse result
    xml_data = results.toxml()
    xml_parsed = BeautifulSoup(xml_data, 'lxml-xml')
    tags = xml_parsed.find_all('literal')
    try:
        structure = tags[0].contents[0]
    except:
        print("No structure found for " + ID)
        wurcs_structure = get_wurcs_from_glytoucan(ID)
        if wurcs_structure:
            structure = glycoct.dumps(wurcs.loads(wurcs_structure))
        else:
            structure = ""
    return structure


if __name__ =='__main__':
    main()
