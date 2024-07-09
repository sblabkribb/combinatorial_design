import pandas as pd
import os
import glob
import pathlib
from dotenv import load_dotenv
load_dotenv()

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Restriction import HincII

from pydna.all import pcr, Assembly, read, Dseqrecord, assembly_fragments
from pydna.primer import Primer
from primers import create
from pydna.common_sub_strings import terminal_overlap

from IPython.display import clear_output

## prepare parts for golden gate assembly ready
def part_insert_goldengate(project_dir, design_file):
    
    """ Generate GenBank files for inserts """
    current_path= pathlib.Path(__file__).parent.absolute()
    project_common_path = os.path.join(current_path, os.getenv("PROJECT_DIR"))
    project_path = os.path.join(project_common_path, project_dir)
    genbank_path = os.path.join(current_path, os.getenv("GENBANK_DIR"))
    part_path = os.path.join(project_path, os.getenv("PART_DIR"))
    insert_path = os.path.join(part_path, "insert")
    design_file_path = os.path.join(project_path, design_file)
    part_dbfile_path = os.path.join(genbank_path, os.getenv("PART_DB_FILE"))


    if not os.path.exists(project_common_path):
        os.makedirs(project_common_path)

    if not os.path.exists(project_path):
        os.makedirs(project_path)

    if not os.path.exists(part_path):
        os.makedirs(part_path)

    ## create the folder for the insert
    if not os.path.exists(insert_path):
        os.makedirs(insert_path)

    # load part design file
    part_data = pd.read_excel(design_file_path)
    # print(part_data.head())
    # exit

    part_data['id'] = part_data['id'].ffill()

    # load part database 
    part_sheets = pd.ExcelFile(part_dbfile_path).sheet_names
    sequences = {sheet: pd.read_excel(part_dbfile_path, sheet_name=sheet) for sheet in part_sheets}

    
    for id in part_data['id'].unique():
        parts = part_data[part_data['id'] == id]
        combined_sequence = ""
        features = []
        # print(id)
        
        # Combine sequences and create features
        start = 0
        for _, part_row in parts.iterrows():
            # print(part_row['name'], part_row['type'])
            part_info = sequences[part_row['type']]
            # print("partinfo:", part_info.head())
            seq = part_info.loc[part_info['Name'] == part_row['name'], 'Sequence'].iloc[0]
            # print("part_info name:", part_info['Name'], "part_row name:", part_row['name'])
            # print(seq)

            if part_row['strand'] == '-':
                seq = str(Seq(seq).reverse_complement())  # Get reverse complement if strand is '-'
            end = start + len(seq)
            combined_sequence += seq
            feature = SeqFeature(FeatureLocation(start=start, end=end), type="misc_feature", qualifiers={"label": part_row['name']})
            features.append(feature)
            start = end

        # Define the GenBank file path
        genbank_path = f'{insert_path}/{id}-insert.gb'

        # Create the SeqRecord
        seq_record = SeqRecord(Seq(combined_sequence), id=id, name=id, description=f"Synthetic construct of {id}")
        seq_record.features = features
        seq_record.annotations = { "molecule_type" : "DNA" }
        
        # print(seq_record)

        # Write to GenBank file
        SeqIO.write(seq_record, genbank_path, 'genbank')

    print(f"Process completed. GenBank files have been created in the folder of {insert_path}")

# Search by label and return the location
def get_positions_by_label(record, labels):
    feature_list = record.features
    startList = []
    endList = []
    for feature in feature_list:
        for label in labels:
            if "label" in feature.qualifiers and feature.qualifiers['label'] == [label]:
                startList.append(int(feature.location.start))
                endList.append(int(feature.location.end))
                # return {'start': int(feature.location.start), 'end': int(feature.location.end)}
    if len(startList) > 0 and len(endList) > 0:
        return {'start': min(startList), 'end': max(endList)}
    else:
        return {'start': None, 'end': None}

## Clone the part insert into a vector using Gibson assembly
def part_withvector_gibson(project_dir, vector_gbfile, target_feature="MCS"):

    current_path= pathlib.Path(__file__).parent.absolute()
    project_common_path = os.path.join(current_path, os.getenv("PROJECT_DIR"))
    project_path = os.path.join(project_common_path, project_dir)
    genbank_path = os.path.join(current_path, os.getenv("GENBANK_DIR"))
    part_path = os.path.join(project_path, os.getenv("PART_DIR"))
    insert_path = os.path.join(part_path, "insert")
    vector_file_path = os.path.join(genbank_path, vector_gbfile)
    withvector_path = os.path.join(part_path, "withvector")

    if not os.path.exists(project_path):
        ## warning no project directory
        print(f"Warning: The project directory '{project_path}' does not exist.")
        return

    ## create the folder for the withvector construct
    if not os.path.exists(f"{withvector_path}"):
        os.makedirs(f"{withvector_path}")

    # read genbank file
    puc19 = read(vector_file_path)

    # no support vector linearization by PCR so cur first and PCR later
    # linearize by enzyme first and then do PCR (simulation only)
    tmp_vector = puc19.linearize(HincII)
    # print(tmp_vector.figure())
    # get primers for vector
    mcs_pos = get_positions_by_label(puc19, [target_feature])
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

    # print(vector_amp.figure())
    # print(vector_amp.list_features())
    primer_list = []

    part_insert_filenames = glob.glob(f"{insert_path}/*-insert.gb")
    for fn in part_insert_filenames:
        insert = Dseqrecord(read(fn))
        # print(insert.seq)
        # print(fn)

        ## primers with homologous region
        # fwd, rev = create(str(insert.seq), add_fwd = Seq(vfwd.seq).reverse_complement(), add_rev = Seq(vrev.seq).reverse_complement())
        ## no homologous region but will be added by Gibson Assembly using assembly_fragments in the following step
        fwd, rev = create(str(insert.seq))

        ifwd = Primer(fwd.seq)
        irev = Primer(rev.seq)
        # print(ifwd)
        # print(irev)
        while len(insert.seq) < len(ifwd) + len(irev):
            print("Warning: The insert sequence is shorter than the primers")
            ifwd = Primer(fwd.seq[:(len(ifwd)-1)])
            irev = Primer(rev.seq[:(len(irev)-1)])
        # print(ifwd)
        # print(irev)

        ifwd.name = "insert_forward"
        irev.name = "insert_reverse"

        insert_amp = pcr(ifwd, irev, insert)
        # print(insert_amp.program())
        # print(insert_amp.figure())
        # print(insert_amp.list_features())
        # print(insert_amp.figure())
        # print(insert_amp.list_features())
        
        vector = Dseqrecord(vector_amp.seq)
        vector.features = vector_amp.features

        # primer update using pcr vector (Dseqrecord type)
        fragment_list = assembly_fragments((vector, insert_amp, vector))
        # print(fragment_list)
        # replace the amplicon with updated primers (have longer homologous region)
        fragment_list = [vector_amp, fragment_list[1]] 

        
        ## primers for vector
        if len(primer_list) == 0:
            primer_list.append({"target": "vector", "direction": "foward", "sequence": str(fragment_list[0].forward_primer.seq), "tm": mt.Tm_NN(fragment_list[0].forward_primer.seq), "length": len(fragment_list[0].forward_primer.seq)})
            primer_list.append({"target": "vector", "direction": "reverse", "sequence": str(fragment_list[0].reverse_primer.seq), "tm": mt.Tm_NN(fragment_list[0].reverse_primer.seq), "length": len(fragment_list[0].reverse_primer.seq)})

        asm = Assembly(fragment_list, algorithm=terminal_overlap)
        candidate = asm.assemble_circular()
        # print(f"total candidates: {len(candidate)}")
        # print(candidate[0].figure())
        # print(candidate[0].list_features())

        ## primers for insert
        primer_list.append({"target": insert.id, "direction": "foward", "sequence": str(fragment_list[1].forward_primer.seq), "tm": mt.Tm_NN(fragment_list[1].forward_primer.seq), "length": len(fragment_list[1].forward_primer.seq)})
        primer_list.append({"target": insert.id, "direction": "reverse", "sequence": str(fragment_list[1].reverse_primer.seq), "tm": mt.Tm_NN(fragment_list[1].reverse_primer.seq), "length": len(fragment_list[1].reverse_primer.seq)})
        
        new_construct = candidate[0]
        # print(type(new_construct))
        

        ## i don't know why, but the feature location is not correct in linux environment
        if os.name == 'posix':
            ## make correction of feature location
            ## this is because the feature shift not applied by pydna
            updated_primer = fragment_list[1].forward_primer.seq
            org_primer = insert_amp.forward_primer.seq
            shift_len = len(updated_primer) - len(org_primer) 
            # print(f"length to shift: {shift_len}")

            insert_pos = get_positions_by_label(new_construct, ["insert_forward"])
            # print(insert_pos)

            # print(new_construct.list_features())
            # new_construct.write("test1.gb")
            # shift all the feature locations from the insert position start
            for feature in new_construct.features:
                if feature.location.start >= insert_pos['start']:
                    # print the direction
                    strand = feature.location.strand
                    # shift the location
                    feature.location = FeatureLocation(start = feature.location.start + shift_len, end = feature.location.end + shift_len)
                    # print the direction
                    feature.location.strand = strand

            # print(new_construct.list_features())         
            # new_construct.write("test2.gb")
    
        # compare all the pairs of features to check the duplication of features. remove all the duplicated features except the first one
        # print(new_construct.list_features())
        feature_duplication = [False] * len(new_construct.features)
        for i in range(len(new_construct.features)):
            for j in range(i+1, len(new_construct.features)):
                if new_construct.features[i].location.start == new_construct.features[j].location.start and new_construct.features[i].location.end == new_construct.features[j].location.end:
                    feature_duplication[j] = True
        new_construct.features = [new_construct.features[i] for i in range(len(new_construct.features)) if not feature_duplication[i]]
        # print(new_construct.list_features())
                
        # print(new_construct.list_features())
        # change the file name to -withvector
        gb_filename = f"{fn.split(os.sep)[-1].replace('-insert.gb', '-withvector.gb')}"
        # gb_filename = f"{fn.split('/')[-1].replace('-insert.gb', '-withvector.gb')}"
        # store the withvector construct
        new_construct.write(f"{withvector_path}/{gb_filename}")
        clear_output(wait=True)
        # clear_output(wait=True)

    # store primer list
    df = pd.DataFrame(primer_list)
    df.to_csv(f"{project_path}/primers-part-withvector-gibson.csv", index=False)
    
    print(f"Process completed. GenBank files have been created in the folder of {withvector_path}")
