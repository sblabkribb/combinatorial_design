import pandas as pd
from itertools import product
import os
import sys
import pathlib
from IPython.display import clear_output
import json

from pydna.common_sub_strings import terminal_overlap
from pydna.common_sub_strings import common_sub_strings
from pydna.all import pcr, Assembly, read, Dseqrecord, assembly_fragments
import part_preparation as pprep
from pydna.primer import Primer
from primers import create
from Bio.Restriction import HincII, BsaI
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqUtils import MeltingTemp as mt

def module_combinatorial_gibson_assembly(project_dir, path_design_file, vector_gbpath, homo="VA"):
    current_path = os.path.dirname(os.path.realpath(__file__))
    pathway_path = os.path.join(current_path, project_dir, "pathway")
    pathway_insert_path = os.path.join(pathway_path, "insert")
    vector_path = os.path.join(current_path, vector_gbpath) 
    # print(pathway_path)

    if not os.path.exists(pathway_path):
        os.makedirs(pathway_path)
        
    if not os.path.exists(pathway_insert_path):
        os.makedirs(pathway_insert_path)

    module_withvector_path = os.path.join(project_dir, "module", "withvector")
    ## read design excel file
    df = pd.read_excel(path_design_file)
    # print(df)
    
    ## genbank file list from module/withvector
    gb_files = [f for f in os.listdir(module_withvector_path) if f.endswith(".gb")]
    # print(gb_files)

    combinatorial_list = {}
    assembly_dictionary = {}
    ## for each row in the design file, find the corresponding genbank files
    for row in df.iterrows():
        ## convert row to series
        row = row[1].dropna()
        print(row[0], "...")

        pathway_mm_path = os.path.join(pathway_insert_path, row[0])

        if not os.path.exists(pathway_mm_path):
            os.makedirs(pathway_mm_path)

        combinations = {}
        ## find the corresponding genbank file that name contains the part name in the gb_files (exclude index column by [1:])
        for r in row[1:]:
            file_list = []
            for gb_file in gb_files:
                if r in gb_file:
                    file_list.append(gb_file)
                
            combinations[r] = file_list
        # print(combinations)

        ## generate all possible combinations
        all_combinations = list(product(*combinations.values()))
        # print(pd.DataFrame(all_combinations))
        
        assembly_list = []
        ## for each combination, do the assembly
        for comb in all_combinations:
            # print(comb)
            combi = all_combinations.index(comb)
            ## print "." for each iteration and new line for every 100 iterations
            ## no new line for the first iteration
            # if combi  % 20 == 0 and combi != 0:
            #     print()
            # print(".", end="")

            ## read the genbank files
            records = [read(os.path.join(module_withvector_path, c)) for c in comb]
            # print(records)

            amp_list = []
            ## for each loaded genbank information, get the VA region for primers and do PCR
            tmp_feature_list = []
            for record in records:
                va_list = []
                for feature in record.features:
                    ## features which labels starting with "VA"
                    # print(feature.qualifiers["label"][0])
                    if feature.qualifiers["label"][0].startswith("VA"):
                        va_list.append(feature.qualifiers["label"][0])
                # print(va_list)
                # print(record.list_features())
                insert_pos = pprep.get_positions_by_label(record, va_list)
                # print(insert_pos)

                frag = record[insert_pos['start']:insert_pos['end']] 
                fwd, rev = create(str(frag.seq),)

                ifwd = Primer(fwd.seq)
                irev = Primer(rev.seq)

                insert_amp = pcr(ifwd, irev, record.seq)
                ## update features
                tmp_feature_list.append(record.features)
                insert_amp.features = frag.features
                amp_list.append(insert_amp)
                # print(frag.list_features())

            ## do the assembly
            asm = Assembly(amp_list, algorithm = terminal_overlap)
            # print(asm)
            candidate = asm.assemble_linear()
            # print(len(candidate))
            # print(candidate[0].figure())
            # print(candidate.list_features())
            ## make 0 lead string for combi
            combistr = str(combi).zfill(4)
            
            candidate[0].name = row[0]+"_"+combistr
            assembly_list.append(candidate[0])

            # write the assembled sequence to a genbank file
            # fname = f"{'_'.join(comb).replace('_withvector.gb', '')}.gb"
            # print(fname)
            combinatorial_list[row[0]+"_"+combistr] = comb
            # print(candidate[0].seq)
            candidate[0].write(os.path.join(pathway_mm_path, candidate[0].name + ".gb"))
            clear_output(wait=True)
        assembly_dictionary[row[0]] = assembly_list

    ## write combinatorial list to a json file
    with open(os.path.join(current_path, project_dir, "module_combinatorial_list.json"), "w") as f:
        json.dump(combinatorial_list, f)
        
    return assembly_dictionary




def module_comb_withvector_gibson_assembly(project_dir, vector_gbpath, linear_insert_dic):

    current_path = os.path.dirname(os.path.realpath(__file__))
    project_path = os.path.join(current_path, project_dir)
    pathway_path = os.path.join(project_path, "pathway")
    pathway_withvector_path = os.path.join(pathway_path, "withvector")
    vector_path = os.path.join(current_path, vector_gbpath) 

    if not os.path.exists(pathway_withvector_path):
        os.makedirs(pathway_withvector_path)
        
    ## ---- assume puc19, need to change if you use other vector ---- ##
    ## read the vector file
    puc19 = read(vector_path)

    # no support vector linearization by PCR so cur first and PCR later
    # linearize by enzyme first and then do PCR (simulation only)
    tmp_vector = puc19.linearize(HincII)
    # print(tmp_vector.figure())
    # get primers for vector
    mcs_pos = pprep.get_positions_by_label(puc19, ["MCS"])
    vector_linear_frag1 = puc19[:mcs_pos['start']]
    vector_linear_frag2 = puc19[mcs_pos['end']:] 
    vfwd1, vrev1 = create(str(vector_linear_frag1.seq))
    vfwd2, vrev2 = create(str(vector_linear_frag2.seq))

    # primers for whole vector 
    vfwd = Primer(vrev1.seq)
    vrev = Primer(vfwd2.seq)
    vfwd.name = "vector_forward"
    vrev.name = "vector_reverse"
    vector_amp = pcr(vfwd, vrev, tmp_vector)

    vector = Dseqrecord(vector_amp.seq)
    vector.features = vector_amp.features
    # print(vector.figure())
    # print(vector.list_features())

    # print(linear_inserts)
    assembly_dictionary = {}
    for module_name, linear_inserts in linear_insert_dic.items():

        pathway_mm_path = os.path.join(pathway_withvector_path, module_name)
        
        if not os.path.exists(pathway_mm_path):
            os.makedirs(pathway_mm_path)

        assembly_list = []    
        primer_list = []
        for insert in linear_inserts:

            amp_list = []
            amp_list.append(vector)

            fwd, rev = create(str(insert.seq),)

            ifwd = Primer(fwd.seq)
            irev = Primer(rev.seq)

            insert_amp = pcr(ifwd, irev, insert.seq)
            ## update features
            # tmp_feature_list.append(record.features)
            insert_amp.features = insert.features
            amp_list.append(insert_amp)
            # print(insert.list_features())
            # print(insert_amp.list_features())

            amp_list.append(vector)

            # print(amp_list)

            ## do the assembly
            fragment_list = assembly_fragments(amp_list)

            # print(fragment_list)
            fragment_list[1].features = insert.features

            ## there is a problem that the features are not shifted properly
            ## the new primer (extended sequence length) of insert_amp is not considered 
            ## split the primer with the terminal 5 base pairs 
            splstr = str(fragment_list[1].primers()[0].seq).split(str(insert.seq)[0:5])
            ## take the first part of the split string and calculate the length
            shift_len = len(splstr[0])
            
            # shift the feature positions to shift_len
            for feature in fragment_list[1].features:
                strand = feature.location.strand
                start = feature.location.start
                end = feature.location.end
                feature.location = FeatureLocation(start = start + shift_len , end = end + shift_len , strand = strand)

            ## add primers for the fragment_list[1]
            primer_list.append({"target": insert.name, "direction": "forward", "sequence": str(fragment_list[1].forward_primer.seq), "tm": mt.Tm_NN(fragment_list[1].forward_primer.seq), "length": len(fragment_list[1].forward_primer.seq)})
            primer_list.append({"target": insert.name, "direction": "reverse", "sequence": str(fragment_list[1].reverse_primer.seq), "tm": mt.Tm_NN(fragment_list[1].reverse_primer.seq), "length": len(fragment_list[1].reverse_primer.seq)})

            # asm = Assembly(amp_list[:-1], algorithm = common_sub_strings, limit=57)
            asm = Assembly(fragment_list)
            # asm = Assembly(amp_list, algorithm = terminal_overlap)
            # print(asm)
            # return 
        
            # candidate = asm.assemble_linear()
            candidate = asm.assemble_circular()
            # print(candidate.synced(vector))
            # print(fragment_list)
            # print(len(candidate))
            # print(candidate[0].figure())
            # print(candidate[0].list_features())

            ## write the assembled sequence to a genbank file
            fname = insert.name + "_withvector.gb"
            # print(fname)

            candidate[0].write(os.path.join(pathway_mm_path, fname))
            clear_output(wait=True)

            assembly_list.append(candidate[0])
        assembly_dictionary[module_name] = assembly_list

    primer_fname = "pathway_primers.csv"
    primer_df = pd.DataFrame(primer_list)
    primer_df.to_csv(os.path.join(project_path, primer_fname), index = False)

    return assembly_dictionary

    ## write the assembled sequence to a genbank file
    # print("writing genbank files...")
    # for i, a in enumerate(assembly_list):
        # fname = "_".join(all_combinations[i].replace("_withvector.gb")) + ".gb"
        # print(fname)
        # a.write(os.path.join(pathway_path, "_".join(all_combinations[i]) + ".gb"))
    


def pathway_comb_withvector_gibson_assembly(project_dir, path_design_file, vector_gbpath, linear_insert_dict):

    current_path = os.path.dirname(os.path.realpath(__file__))
    project_path = os.path.join(current_path, project_dir)
    pathway_path = os.path.join(project_path, "pathway")
    pathway_withvector_path = os.path.join(pathway_path, "withvector")
    vector_path = os.path.join(current_path, vector_gbpath) 

    if not os.path.exists(pathway_withvector_path):
        os.makedirs(pathway_withvector_path)
        
    ## read design excel file from the second sheet
    df = pd.read_excel(path_design_file, sheet_name=1)

    ## for each row in the design file, find the corresponding genbank information from linear_insert_dict
    for row in df.iterrows():
        ## convert row to series
        row = row[1].dropna()
        module_id = row[0]
        print(module_id, "...")
        
        pathway_mm_path = os.path.join(pathway_withvector_path, module_id)

        if not os.path.exists(pathway_mm_path):
            os.makedirs(pathway_mm_path)

        ## filtering out linear_insert_dict items which are not specified in the df columns 
        linear_insert_dict = {k: v for k, v in linear_insert_dict.items() if k in row[1:].values}

    ## ---- assume puc19, need to change if you use other vector ---- ##
    ## read the vector file
    puc19 = read(vector_path)

    # no support vector linearization by PCR so cur first and PCR later
    # linearize by enzyme first and then do PCR (simulation only)
    tmp_vector = puc19.linearize(HincII)
    # print(tmp_vector.figure())
    # get primers for vector
    mcs_pos = pprep.get_positions_by_label(puc19, ["MCS"])
    vector_linear_frag1 = puc19[:mcs_pos['start']]
    vector_linear_frag2 = puc19[mcs_pos['end']:] 
    vfwd1, vrev1 = create(str(vector_linear_frag1.seq))
    vfwd2, vrev2 = create(str(vector_linear_frag2.seq))

    # primers for whole vector 
    vfwd = Primer(vrev1.seq)
    vrev = Primer(vfwd2.seq)
    vfwd.name = "vector_forward"
    vrev.name = "vector_reverse"
    vector_amp = pcr(vfwd, vrev, tmp_vector)

    vector = Dseqrecord(vector_amp.seq)
    vector.features = vector_amp.features
    # print(vector.figure())
    # print(vector.list_features())

    ## find the corresponding keys from the linear_insert_dict 

    ## all combinations of linear inserts from the linear_insert_dict
    all_combinations = list(product(*linear_insert_dict.values()))

    assembly_list_insertonly = []
    ## assembly the pathways
    for i, linear_inserts in enumerate(all_combinations):

        ## print "." for each iteration and new line for every 100 iterations
        ## no new line for the first iteration
        if i  % 50 == 0 and i != 0:
            print(i)
            print()
        print(".", end="")
        
        ## assemble the modules into a linear dna 
        amp_list = []
        for insert in linear_inserts:
            fwd, rev = create(str(insert.seq))

            ifwd = Primer(fwd.seq)
            irev = Primer(rev.seq)

            insert_amp = pcr(ifwd, irev, insert.seq)
            ## update features
            insert_amp.features = insert.features
            amp_list.append(insert_amp)
            # print(insert.list_features())

        asm = Assembly(amp_list, algorithm = terminal_overlap)
        # print(asm)
        candidate = asm.assemble_linear()
        # print(len(candidate))
        # print(candidate[0].figure())
        # print(candidate.list_features())
        ## make 0 lead string for combi
        # combistr = str(combi).zfill(4)
        
        candidate[0].name = module_id + "_" + str(i).zfill(4)
        assembly_list_insertonly.append(candidate[0])
        if i == 10:
            break

    ## assembly with the insert and vector
    assembly_list_final = []
    primer_list = []
    for insert in assembly_list_insertonly:

        fwd, rev = create(str(insert.seq))

        ifwd = Primer(fwd.seq)
        irev = Primer(rev.seq)

        insert_amp = pcr(ifwd, irev, insert.seq)
        ## update features
        # tmp_feature_list.append(record.features)
        insert_amp.features = insert.features
        # print(amp_list)

        ## do the assembly
        frag_list = [vector, insert_amp, vector]
        fragment_list = assembly_fragments(frag_list)

        # print(fragment_list)
        fragment_list[1].features = insert.features

        ## there is a problem that the features are not shifted properly
        ## the new primer (extended sequence length) of insert_amp is not considered 
        ## split the primer with the terminal 5 base pairs 
        splstr = str(fragment_list[1].primers()[0].seq).split(str(insert.seq)[0:5])
        ## take the first part of the split string and calculate the length
        shift_len = len(splstr[0])
        
        # shift the feature positions to shift_len
        for feature in fragment_list[1].features:
            strand = feature.location.strand
            start = feature.location.start
            end = feature.location.end
            feature.location = FeatureLocation(start = start + shift_len , end = end + shift_len , strand = strand)

        ## add primers for the fragment_list[1]
        primer_list.append({"target": insert.name, "direction": "forward", "sequence": str(fragment_list[1].forward_primer.seq), "tm": mt.Tm_NN(fragment_list[1].forward_primer.seq), "length": len(fragment_list[1].forward_primer.seq)})
        primer_list.append({"target": insert.name, "direction": "reverse", "sequence": str(fragment_list[1].reverse_primer.seq), "tm": mt.Tm_NN(fragment_list[1].reverse_primer.seq), "length": len(fragment_list[1].reverse_primer.seq)})

        asm = Assembly(fragment_list)
        # asm = Assembly(amp_list, algorithm = terminal_overlap)
        # candidate = asm.assemble_linear()
        candidate = asm.assemble_circular()
        # print(candidate.synced(vector))
        # print(fragment_list)
        # print(len(candidate))
        # print(candidate[0].figure())
        # print(candidate[0].list_features())

        # a.write(os.path.join(pathway_mm_path, a.name + ".gb"))
        # clear_output(wait=True)
        # print(amp_list)

        ## write the assembled sequence to a genbank file
        fname = insert.name + "_withvector.gb"
        # print(fname)

        candidate[0].write(os.path.join(pathway_mm_path, fname))
        clear_output(wait=True)

        assembly_list_final.append(candidate[0])

    primer_fname = "multi_pathway_primers.csv"
    primer_df = pd.DataFrame(primer_list)
    primer_df.to_csv(os.path.join(project_path, primer_fname), index = False)

    return assembly_list_final

    ## write the assembled sequence to a genbank file
    # print("writing genbank files...")
    # for i, a in enumerate(assembly_list):
        # fname = "_".join(all_combinations[i].replace("_withvector.gb")) + ".gb"
        # print(fname)
        # a.write(os.path.join(pathway_path, "_".join(all_combinations[i]) + ".gb"))


