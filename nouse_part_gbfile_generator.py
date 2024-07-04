
curpath = "/home/haseong/dev/mvaopt/labnote/002_part_preparation"

# read genbank file
puc19 = read(f"{curpath}/genbank/addgene-plasmid-50005-sequence-222046.gb")

# no support vector linearization by PCR so cur first and PCR later
# linearize by enzyme first and then do PCR (simulation only)
tmp_vector = puc19.linearize(HincII)
# print(tmp_vector.figure())
# get primers for vector
mcs_pos = get_positions_by_label(puc19, "MCS")
vector_linear_frag1 = puc19[:mcs_pos['start']]
vector_linear_frag2 = puc19[mcs_pos['end']:] 
vfwd1, vrev1 = create(str(vector_linear_frag1.seq))
vfwd2, vrev2 = create(str(vector_linear_frag2.seq))

# PCR whole vector with proper primers 
vfwd = Primer(vrev1.seq)
vrev = Primer(vfwd2.seq)
vfwd.name = "vector_forward"
vrev.name = "vector_reverse"
vector_amp = pcr(vfwd, vrev, tmp_vector)

# print(vector_amp.figure())
# print(vector_amp.list_features())

primer_list = []

part_insert_filenames = glob.glob(f"{curpath}/genbank_part/*-insert.gb")
for fn in part_insert_filenames:
    insert = Dseqrecord(read(fn))
    # print(insert.seq)
    print(fn)

    ## primers with homologous region
    # fwd, rev = create(str(insert.seq), add_fwd = Seq(vfwd.seq).reverse_complement(), add_rev = Seq(vrev.seq).reverse_complement())
    ## no homologous region but will be added by Gibson Assembly using assembly_fragments in the following step
    fwd, rev = create(str(insert.seq))

    ifwd = Primer(fwd.seq)
    irev = Primer(rev.seq)
    ifwd.name = "insert_forward"
    irev.name = "insert_reverse"

    insert_amp = pcr(ifwd, irev, insert)
    # print(insert_amp.program())
    # print(insert_amp.figure())
    # print(insert_amp.list_features())

    # print(type(insert))

    ## primers with homologous region
    # fwd, rev = create(str(insert.seq), add_fwd = Seq(vfwd).reverse_complement(), add_rev = Seq(vrev).reverse_complement())
    # tmpdic = {"target": insert.id, "target seq": str(insert.seq), "direction": "foward", "sequence": str(fwd.seq), "tm": fwd.tm, "length": len(fwd.seq)}
    # primer_list.append(tmpdic)
    # tmpdic = {"target": insert.id, "target seq": str(insert.seq), "direction": "reverse", "sequence": str(rev.seq), "tm": rev.tm, "length": len(rev.seq)}
    # primer_list.append(tmpdic)

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
        primer_list.append({"target": "vector", "target seq": str(vector.seq), "direction": "foward", "sequence": str(fragment_list[0].forward_primer.seq), "tm": mt.Tm_NN(fragment_list[0].forward_primer.seq), "length": len(fragment_list[0].forward_primer.seq)})
        primer_list.append({"target": "vector", "target seq": str(vector.seq), "direction": "reverse", "sequence": str(fragment_list[0].reverse_primer.seq), "tm": mt.Tm_NN(fragment_list[0].reverse_primer.seq), "length": len(fragment_list[0].reverse_primer.seq)})


    # this is the way suggesting by the author
    # fragment_list = assembly_fragments((vector, insert_amp, vector))
    # fragment_list = fragment_list[:-1]
    # fragment_list = [vector, insert_amp]
    # asm = Assembly(fragment_list)
    # print(fragment_list[0].figure())
    # print(fragment_list[1].figure())

    asm = Assembly(fragment_list, algorithm=terminal_overlap)
    candidate = asm.assemble_circular()
    print(f"total candidates: {len(candidate)}")
    # print(candidate[0].figure())
    # print(candidate[0].list_features())

    ## primers for insert
    primer_list.append({"target": insert.id, "target seq": str(insert.seq), "direction": "foward", "sequence": str(fragment_list[1].forward_primer.seq), "tm": mt.Tm_NN(fragment_list[1].forward_primer.seq), "length": len(fragment_list[1].forward_primer.seq)})
    primer_list.append({"target": insert.id, "target seq": str(insert.seq), "direction": "reverse", "sequence": str(fragment_list[1].reverse_primer.seq), "tm": mt.Tm_NN(fragment_list[1].reverse_primer.seq), "length": len(fragment_list[1].reverse_primer.seq)})
    
    new_construct = candidate[0]

    ## make a correction because the feature shift not applied by pydna
    updated_primer = fragment_list[1].forward_primer.seq
    org_primer = insert_amp.forward_primer.seq
    shift_len = len(updated_primer) - len(org_primer) 
    # print(f"length to shift: {shift_len}")

    insert_pos = get_positions_by_label(new_construct, "insert_forward")
    # print(insert_pos)

    # # shift all the feature locations from the insert position start
    for feature in new_construct.features:
        if feature.location.start >= insert_pos['start']:
            # print the direction
            strand = feature.location.strand
            # shift the location
            feature.location = FeatureLocation(start = feature.location.start + shift_len, end = feature.location.end + shift_len)
            # print the direction
            feature.location.strand = strand
            

    # print(new_construct.list_features())
    gb_filename = fn.replace("-insert", "-withvector")
    # print(gb_filename)
    # break
    new_construct.write(gb_filename)

# store primer list
df = DataFrame(primer_list)
df.to_csv(f"{curpath}/primers.csv", index=False)
    
