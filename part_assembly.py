import pandas as pd
from itertools import product
import os
import pathlib

from pydna.common_sub_strings import terminal_overlap
from pydna.all import pcr, Assembly, read, Dseqrecord, assembly_fragments
import part_preparation as pprep
from pydna.primer import Primer
from primers import create
from Bio.Restriction import HincII, BsaI
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqUtils import MeltingTemp as mt

from IPython.display import clear_output


def get_part_combinations(design_file):

    ## current path
    current_path = os.path.dirname(os.path.realpath(__file__))

    assembly_design_path = f"{current_path}/{design_file}"
    # Load the Excel file
    df = pd.read_excel(assembly_design_path)
    df['id'] = df['id'].ffill()

    # for those that have the same id
    module_ids = df['id'].unique()
    # print(module_ids)

    combinations = []

    # read part type for each id
    for id in module_ids:
        # get the parts for the module
        parts = df[df['id'] == id]
        part_types = parts['type']
        # print(part_types)
        
        # Dictionary to hold the parts by type
        parts_by_type = {}
        
        new_part_types = [] 
        for i,part_type in enumerate(part_types):
        # for part_type in part_types:
            # Collect all non-null items for each type across all relevant columns
            # tries to 
            new_type = part_type+"_"+str(i)
            new_part_types.append(new_type)
            parts_by_type[new_type] = parts.iloc[i]["name":].values.flatten()
            # print(parts.loc[parts['type'] == part_type, 'name':].values.flatten())
            # parts_by_type[part_type] = parts.loc[parts['type'] == part_type, 'name':].values.flatten()
            parts_by_type[new_type] = [part for part in parts_by_type[new_type] if pd.notna(part)]

        # print(parts_by_type)

        # Generate all combinations of the specified part types
        all_combinations = list(product(*[parts_by_type[new_type] for new_type in new_part_types]))

        # Convert each combination into a dictionary
        dictionaries = []
        for combination in all_combinations:
            dictionary = {part_type: combination[i] for i, part_type in enumerate(new_part_types)}
            dictionaries.append(dictionary)

        combinations.extend(dictionaries)

    return combinations


## assembly the parts from all the combinations of the list
def get_all_possible_combinations(project_dir, combinations):
    
    ## list all the parts in the withvector directory
    current_path = os.path.dirname(os.path.realpath(__file__))
    withvector_path = os.path.join(current_path, project_dir, "part", "withvector")
    withvector_files = os.listdir(withvector_path)
    # print(withvector_files)
    
    # a dictionary of all the combinations
    all_combinations = {}

    ## for each combination, find the parts in the withvector directory
    ## combinations is a list of dictionaries
    for combination in combinations:
        # print(combination)
        # loop through the values and keys at the same time 
        # a dictionary for storing the files for each part
        combination_dict = {}
        for part_type, part_name in combination.items():
            # search for all the files where the part_name is in the file name
            # a list storing the files for each part
            part_list = []
            for withvector_file in withvector_files:
                if part_name in withvector_file:
                    part_list.append(withvector_file)
                    # print(f"Found part {part_name} in withvector directory")
                    # print(withvector_file)
                    # break
            combination_dict[part_type] = part_list
        # print(combination_dict)

        # generate a index key for the combination
        index_key = "-".join([part_name for part_name in combination.values()])
        # print(index_key)
        # the length of index_key should be the same as the number of part combinations specified in the excel file
         
        # Generate all combinations of the specified gb files
        all_combinations[index_key] = (list(product(*[combination_dict[part_type] for part_type in combination_dict])))


    return all_combinations


def filtering_combinations(all_combination_list):
    
    ## all possible combinations with a two column list
    # transform the dictionary into a list with two columns "key" and "value"
    # if value is a list, then expand the list
    all_combinations = []
    for key, value in all_combination_list.items():
        if isinstance(value, list):
            for item in value:
                all_combinations.append([key, item])
        else:
            all_combinations.append([key, value])
    
    ## filtering out if its overhang junction is not the same between two parts
    filteredList = {}

    for i in range(len(all_combinations)): 
        ## seperate the list into two groups where one is 0 to len-1 and the other is 1 to len
        group1 = all_combinations[i][1][:-1]
        group2 = all_combinations[i][1][1:]
        # print(group1)
        ## separate each items with "-" in the groups
        overhang1 = [s.replace("-withvector.gb", "").split("-")[-1] for s in group1]
        overhang2 = [s.replace("-withvector.gb", "").split("-")[0] for s in group2]
        # print(overhang1)
        # print(overhang2)
        ## then compare the overhang junctions between two parts
        if(overhang1 == overhang2):
            filteredList[all_combinations[i][0]]=all_combinations[i][1]
            # go for the assembly

    return filteredList


def part_assembly_goldengate(project_dir, filtered_combinations):
    # current directory
    current_dir = pathlib.Path(__file__).parent.absolute()
    withvector_path = os.path.join(current_dir, project_dir, "part", "withvector") 
    module_dir = os.path.join(current_dir, project_dir, "module")
    module_insert_dir = os.path.join(module_dir, "inserts")

    if not os.path.exists(module_dir):
        os.makedirs(module_dir)

    if not os.path.exists(module_insert_dir):
        os.makedirs(module_insert_dir)
    
    candidates = []
    # read the gb files in filtered_combinations from withvector directory
    for key, value in filtered_combinations.items():
        # print(value)
        fragments = []
        for gbfile in value:
            record = read(f"{withvector_path}/{gbfile}")
            # print(record.list_features())
            # print(gbfile)

            part_of_interest = gbfile.replace("-withvector.gb", "").split("-")
            # POI plus BsaI site
            part_of_interest.append("BsaI")
            pos = pprep.get_positions_by_label(record, part_of_interest)
            # print(pos)
            record_sub = record[pos['start']:pos['end']]
            record_seq = record_sub.seq

            ## --- feature location ---##
            # list all the features where start is greater than pos['start'] and end is less than pos['end']
            record_sub.features = [feature for feature in record.features if feature.location.start >= pos['start'] and feature.location.end <= pos['end']]
            # print(record_sub.list_features())

            ## find features type is primer_bind
            feature_index = []
            for feature in record_sub.features:
                if feature.type == "primer_bind":
                    # print(feature.qualifiers["label"])
                    feature_index.append(record_sub.features.index(feature))
            
            ## remove features from the reverse index 
            ## this is because the features are not removed properly if we remove from the first index
            for i in feature_index[::-1]:
                record_sub.features.pop(i)
            # print(record_sub.list_features())
            ## ---------------------- ##

            # generate primers for the part
            fwd, rev = create(str(record_seq))
            # print(fwd, rev)

            ifwd = Primer(fwd.seq)
            irev = Primer(rev.seq)
            # print(ifwd)
            # print(irev)
            ifwd.name = "insert_forward"
            irev.name = "insert_reverse"

            insert_amp = pcr(ifwd, irev, record_seq)
            insert_amp.features = record_sub.features
            # print(insert_amp)
            # print(insert_amp.list_features())
            
            # find the minimum position of the features
            min_pos1 = min([feature.location.start for feature in insert_amp.features])
            min_pos2 = min([feature.location.end for feature in insert_amp.features])
            min_pos = min(min_pos1, min_pos2)

            # shift the feature positions to the min_pos 
            for feature in insert_amp.features:
                strand = feature.location.strand
                start = feature.location.start
                end = feature.location.end
                feature.location = FeatureLocation(start = start - min_pos , end = end - min_pos , strand = strand)

            # print(insert_amp.list_features())
            # print(insert_amp.program())
            insert_amp_cut = insert_amp.cut(BsaI)
            # print("cut", insert_amp_cut)
            
            # find the longest fragment in the insert_amp_cut
            insert_amp_cut = max(insert_amp_cut, key = lambda x: len(x))

            fragments.append(insert_amp_cut)
            # print(insert_amp_cut.figure())
        
        # print(fragments)
        # print(fragments[0].list_features())
        # print(fragment_list[1])
        # print(fragment_list[1].figure())
        asm = Assembly(fragments, algorithm=terminal_overlap, limit = 4)
        candidate = asm.assemble_linear()
        # print(len(candidate))
        # print(candidate[0].figure())
        # print(candidate[0].list_features())

        ## write the candidate to the module "insert" directory
        fname = f"{module_insert_dir}/{key}.gb"
        candidate[0].write(fname)
        candidate[0].name = fname
        candidates.append(candidate)
        clear_output(wait=True)

    return candidates

    

def build_module_withvector_gibson(project_dir, vector_gbpath, module_linear):
    # current directory
    current_dir = pathlib.Path(__file__).parent.absolute()
    vector_path = os.path.join(current_dir, vector_gbpath) 
    module_dir = os.path.join(current_dir, project_dir, "module")
    module_insert_dir = os.path.join(module_dir, "inserts")
    module_withvector_dir = os.path.join(module_dir, "withvector")

    if not os.path.exists(module_withvector_dir):
        os.makedirs(module_withvector_dir)
    
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

    ## read all the genbank files in the module_insert_dir
    module_insert_files = os.listdir(module_insert_dir)
    # print(module_insert_files)

    candidates = []
    
    primer_list = []

    # read name of each of module_linear elements and find the corresponding file in the module_insert_dir
    for module in module_linear:
        module_linear_file = module[0].name
        # extract file name (remove path)
        
        # print(module_linear_file)
        record = read(module_linear_file)
        # print(record.list_features())
        fwd, rev = create(str(record.seq))
        # print(fwd)

        ifwd = Primer(fwd.seq)
        irev = Primer(rev.seq)
        # ifwd.name = "insert_forward"
        # irev.name = "insert_reverse"

        insert_amp = pcr(ifwd, irev, record.seq)
        insert_amp.features = record.features
        # print(insert_amp)
        # print(insert_amp.list_features())

        # primer update using pcr vector (Dseqrecord type)
        fragment_list = assembly_fragments((vector, insert_amp, vector))
        # replace the amplicon with updated primers (have longer homologous region)
        fragment_list[1].features = insert_amp.features
        # print(fragment_list[1].primers()[0].seq)

        ## there is a problem that the features are not shifted properly
        ## the new primer (extended sequence length) of insert_amp is not considered 
        ## split the primer with the terminal 5 base pairs 
        splstr = str(fragment_list[1].primers()[0].seq).split(str(record.seq)[0:5])
        ## take the first part of the split string and calculate the length
        shift_len = len(splstr[0])
        
        # shift the feature positions to shift_len
        for feature in fragment_list[1].features:
            strand = feature.location.strand
            start = feature.location.start
            end = feature.location.end
            feature.location = FeatureLocation(start = start + shift_len , end = end + shift_len , strand = strand)

        # aln = pairwise2.align.globalxx(str(fragment_list[1].primers()[0].seq), str(record.seq))

        ## add primers for the fragment_list[1]
        primer_list.append({"target": os.path.basename(module_linear_file), "direction": "forward", "sequence": str(fragment_list[1].forward_primer.seq), "tm": mt.Tm_NN(fragment_list[1].forward_primer.seq), "length": len(fragment_list[1].forward_primer.seq)})
        primer_list.append({"target": os.path.basename(module_linear_file), "direction": "reverse", "sequence": str(fragment_list[1].reverse_primer.seq), "tm": mt.Tm_NN(fragment_list[1].reverse_primer.seq), "length": len(fragment_list[1].reverse_primer.seq)})
        
        # print(insert_amp_cut)
        asm = Assembly(fragment_list, algorithm=terminal_overlap)
        # print(asm)
        candidate = asm.assemble_circular()
        # print(len(candidate))
        # print(candidate[0].figure())
        # print(candidate[0].list_features())

        fname = f"{module_withvector_dir}/{os.path.basename(module_linear_file).replace('.gb', '')}_withvector.gb"
        candidate[0].write(fname)
        clear_output(wait=True)

        candidates.append(candidate[0])

    primer_fname = f"{module_dir}/module_primers.csv"
    primer_df = pd.DataFrame(primer_list)
    primer_df.to_csv(primer_fname, index = False)

    return candidates