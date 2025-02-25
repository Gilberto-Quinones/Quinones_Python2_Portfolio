# Quinones_Python2_Portfolio
This is the Portfolio for python2 winter quarter 
Felt like Python Two was more neat and cool 
Open cv was pretty cool 



## This is the Sequence Objects pts 1-4 
```
      from Bio.Seq import Seq
```
      my_seq = Seq("GATCG")
``` 
      for index, letter in enumerate(my_seq)
``` 
      print("%i %s" % (index,letter))
````
      0 G
      1 A
      2 T
      3 C
      4 G
````
# we can also print the length of each sequence 
```
        print(len(my_seq))
```
        5
      
        print(my_seq[0])
```  
        G
```
        print(my_seq[4])
```
      
        G
```
        print(my_seq[2])
```
        T
```
        Seq("AAAA").count("AA")
```
        2
```
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
len(my_seq)
32
my_seq.count("G")
9
100 * (my_seq.count("G") + my_seq.count("C")) / len(my_seq)
46.875
from Bio.SeqUtils import gc_fraction
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
gc_fraction(my_seq)
0.46875
my_seq[4:12]
Seq('GATGGGCC')
my_seq[0::3]
Seq('GCTGTAGTAAG')
my_seq[1::3]
Seq('AGGCATGCATC')
my_seq[2:3]
Seq('T')
my_seq[::-1]
Seq('CGCTAAAAGCTAGGATATATCCGGGTAGCTAG')
str(my_seq)
'GATCGATGGGCCTATATAGGATCGAAAATCGC'
fasta_format_string = ">Name\n%s\n" % my_seq
print(fasta_format_string)
>Name
GATCGATGGGCCTATATAGGATCGAAAATCGC

seq1 = Seq("ACGT")
seq2 = Seq("AACCGG")
seq1 + seq2
Seq('ACGTAACCGG')
seq2 + seq1 
Seq('AACCGGACGT')
contigs = [Seq("ATG"), Seq("ATCCCG"), Seq("TTGA")]
spacer = Seq("N" *10)
spacer.join(contigs)
Seq('ATGNNNNNNNNNNATCCCGNNNNNNNNNNTTGA')
dna_seq = Seq("acgtACGT")
dna_seq
Seq('acgtACGT')
dna_seq.upper()
Seq('ACGTACGT')
dna_seq.lower()
Seq('acgtacgt')
dna_seq.upper()
Seq('ACGTACGT')
"gtac" in dna_seq
False
"GTAC" in dna_seq
False
dna_seq = dna_seq.upper()
"GTAC" in dna_seq
True
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
my_seq.complement()
Seq('CTAGCTACCCGGATATATCCTAGCTTTTAGCG')
my_seq.reverse_complement()
Seq('GCGATTTTCGATCCTATATAGGCCCATCGATC')
protein_seq = Seq("EVRNAK")
protein_seq.complement()
Seq('EBYNTM')
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
coding_dna
Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')
template_dna = coding_dna.reverse_complement()
template_dna
Seq('CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT')
coding_dna
Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')
messenger_rna = coding_dna.transcribe()
messenger_rna
Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')
template_dna.reverse_complement().transcribe()
Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')
messenger_rna.back_transcribe()
Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')
messenger_rna.translate()
Seq('MAIVMGR*KGAR*')
coding_dna.translate(table="Vertebrate Mitochondrial")
Seq('MAIVMGRWKGAR*')
coding_dna.translate(table = 2)
Seq('MAIVMGRWKGAR*')
coding_dna.translate(to_stop = True)
Seq('MAIVMGR')
coding_dna.translate(table = 2, to_stop=True)
Seq('MAIVMGRWKGAR')
coding_dna.translate(table = 2, stop_symbol = "!")
Seq('MAIVMGRWKGAR!')
gene = Seq("ATGGTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTTCTCCCATGGCAGCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGATAATCGTGGCTATTACTGGGUGATAA")
gene.translate(table = "Bacterial")
Seq('MVKKMQSIVLALSLVLVLPWQHRLRKLR*SRQ*NYR*AIVIIVAITG**')
gene.translate(table = "Bacterial", to_stop = True)
Seq('MVKKMQSIVLALSLVLVLPWQHRLRKLR')
gene.translate(table = "Bacterial" , cds = True)




## Sequence Annotation 

```
        from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
simple_seq = Seq("GATC")
simple_seq_r = SeqRecord(simple_seq)
simple_seq_r
SeqRecord(seq=Seq('GATC'), id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=[])
simple_seq_r.id = "AC12345"
simple_seq_r.description = "Made up sequence for the VDB Computational Biology Class"
print(simple_seq_r.description)
Made up sequence for the VDB Computational Biology Class
simple_seq_r.seq
Seq('GATC')
simple_seq_r
SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequence for the VDB Computational Biology Class', dbxrefs=[])
simple_seq_r.annotations["evidence"] ="None. This is just an example "
print(simple_seq_r.annotations["evidence"]) 
None. This is just an example 
simple_seq_r
SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequence for the VDB Computational Biology Class', dbxrefs=[])
simple_seq_r.letter_annotations["phred_quality"]= [40, 40, 39, 30]
print(simple_seq_r.letter_annotations)
{'phred_quality': [40, 40, 39, 30]}
#https://raw.githubbusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.fna
from Bio import SeqIO
record = SeqIO.read("NC_005816.fna.txt", "fasta")
record
SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='gi|45478711|ref|NC_005816.1|', name='gi|45478711|ref|NC_005816.1|', description='gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])
record.seq
Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')
record.id
'gi|45478711|ref|NC_005816.1|'
record.description
'gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'
record.dbxrefs
[]
record.annotations
{}
record.features
[]
#http:raw.githubusercontent.com/biopython/biopython/master/Test/GenBank/NC_005816.gb
record = SeqIO.read("NC_005816.gb.txt", "genbank")
record
SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])
record.seq
Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')
record.id
'NC_005816.1'
record.name
'NC_005816'
record.description
'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'
record.letter_annotations
{}
len(record.annotations)
13
record.annotations["source"]
'Yersinia pestis biovar Microtus str. 91001'
record.dbxrefs
['Project:58037']
len(record.features)
41
from Bio import SeqFeature
start_pos = SeqFeature.AfterPosition(5)
end_pos = SeqFeature.BetweenPosition(9, left=8, right=9)
my_location = SeqFeature.SimpleLocation(start_pos, end_pos)
print(my_location)
[>5:(8^9)]
my_location.end
BetweenPosition(9, left=8, right=9)
int(my_location.end)
9
int(my_location.start)
5
exact_location = SeqFeature.SimpleLocation(5,9)
print(exact_location)
[5:9]
exact_location.start
ExactPosition(5)
from Bio.SeqRecord import SeqRecord
record = SeqRecord(Seq("MMYQQGCFAGGTVRLAKDLAENNRGARVLVVCSEITAVTFRGSPSETHLDSMVGQALFGD"
                "GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK"
                "NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM"
                "SSAC"),
                        id= "gi|14150838|gb |AAK54648.1| AF376133_1",
            description= "chacone synthase [Cucumis sativus]",
                  )
print(record.format("fasta"))
>gi|14150838|gb |AAK54648.1| AF376133_1 chacone synthase [Cucumis sativus]
MMYQQGCFAGGTVRLAKDLAENNRGARVLVVCSEITAVTFRGSPSETHLDSMVGQALFGD
GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK
NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM
SSAC

print(record)
ID: gi|14150838|gb |AAK54648.1| AF376133_1
Name: <unknown name>
Description: chacone synthase [Cucumis sativus]
Number of features: 0
Seq('MMYQQGCFAGGTVRLAKDLAENNRGARVLVVCSEITAVTFRGSPSETHLDSMVG...SAC')
from Bio import SeqIO
record = SeqIO.read("NC_005816.gb.txt", "genbank")
record
SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])
len(record)
9609
len(record.features)
41
print(record.features[20])
type: gene
location: [4342:4780](+)
qualifiers:
    Key: db_xref, Value: ['GeneID:2767712']
    Key: gene, Value: ['pim']
    Key: locus_tag, Value: ['YP_pPCP05']

print(record.features[21])
type: CDS
location: [4342:4780](+)
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
    Key: gene, Value: ['pim']
    Key: locus_tag, Value: ['YP_pPCP05']
    Key: note, Value: ['similar to many previously sequenced pesticin immunity protein entries of Yersinia pestis plasmid pPCP, e.g. gi| 16082683|,ref|NP_395230.1| (NC_003132) , gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655| emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1| (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)']
    Key: product, Value: ['pesticin immunity protein']
    Key: protein_id, Value: ['NP_995571.1']
    Key: transl_table, Value: ['11']
    Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']

sub_record = record[4300:4800]
len(sub_record)
500
len(sub_record.features)
2
sub_record.features[0]
SeqFeature(SimpleLocation(ExactPosition(42), ExactPosition(480), strand=1), type='gene', qualifiers=...)
sub_record.features[1]
SeqFeature(SimpleLocation(ExactPosition(42), ExactPosition(480), strand=1), type='CDS', qualifiers=...)
print(sub_record.features[0])
type: gene
location: [42:480](+)
qualifiers:
    Key: db_xref, Value: ['GeneID:2767712']
    Key: gene, Value: ['pim']
    Key: locus_tag, Value: ['YP_pPCP05']

print(sub_record.features[1])
type: CDS
location: [42:480](+)
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
    Key: gene, Value: ['pim']
    Key: locus_tag, Value: ['YP_pPCP05']
    Key: note, Value: ['similar to many previously sequenced pesticin immunity protein entries of Yersinia pestis plasmid pPCP, e.g. gi| 16082683|,ref|NP_395230.1| (NC_003132) , gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655| emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1| (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)']
    Key: product, Value: ['pesticin immunity protein']
    Key: protein_id, Value: ['NP_995571.1']
    Key: transl_table, Value: ['11']
    Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']

sub_record.annotations
{'molecule_type': 'DNA'}
sub_record.dbxrefs
[]
sub_record.annotations["topology"] = "linear"
sub_record.annotations
{'molecule_type': 'DNA', 'topology': 'linear'}
sub_record.id
'NC_005816.1'
sub_record.name
'NC_005816'
sub_record.description
'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'
sub_record.description = 'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial sequence'
​
sub_record.description
'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial sequence'
print(sub_record.format("genbank")[:200]+ "...")
LOCUS       NC_005816                500 bp    DNA     linear   UNK 01-JAN-1980
DEFINITION  Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial
            sequence.
ACCESSION   NC_00581...
record = SeqIO.read("NC_005816.gb.txt", "genbank")
record
SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])
len(record)
9609
len(record.features)
41
record.dbxrefs
['Project:58037']
record.annotations.keys()
dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])
shifted = record[2000:] + record[:2000]
shifted
SeqRecord(seq=Seq('GATACGCAGTCATATTTTTTACACAATTCTCTAATCCCGACAAGGTCGTAGGTC...GGA'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])
len(shifted)
9609
len(shifted.features)
40
shifted.annotations.keys()
dict_keys(['molecule_type'])
shifted.dbxrefs
[]
shifted.dbxrefs = record.dbxrefs[:]
shifted.dbxrefs
['Project:58037']
shifted.annotations = record.annotations.copy()
shifted.annotations.keys()
dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])
record
SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])
print("%s %i %i %i %i" % (record.id,len(record), len(record.features), len(record.dbxrefs), len(record.annotations)))
NC_005816.1 9609 41 1 13
rc = record.reverse_complement(id = "Testing")
rc
SeqRecord(seq=Seq('CAGGGGTCGGGGTACGCATTCCCTCATGCGTCAATATTATCTGGCATTGCGATG...ACA'), id='Testing', name='<unknown name>', description='<unknown description>', dbxrefs=[])
print("%s %i %i %i %i" % (rc.id, len(rc), len(rc.features), len(rc.dbxrefs), len(rc.annotations)))
Testing 9609 41 0 0

```

# SequenceI/O Section 

```

    .com/biopython/biobython/master/Doc/examples/PF05371_seed.sth
from Bio import AlignIO
alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
print(alignment)
Alignment with 7 rows and 52 columns
AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRL...SKA COATB_BPIKE/30-81
AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKL...SRA Q9T0Q8_BPIKE/1-52
DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRL...SKA COATB_BPI22/32-83
AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73
alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
print("Alignment length %i" %alignment.get_alignment_length())
Alignment length 52
for record in alignment: 
    print("%s - %s" % (record.seq, record.id))
AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73
for record in alignment:
    if record.dbxrefs:
        print("%s %s" % (record.id, record.dbxrefs))
COATB_BPIKE/30-81 ['PDB; 1ifl ; 1-52;']
COATB_BPM13/24-72 ['PDB; 2cpb ; 1-49;', 'PDB; 2cps ; 1-49;']
Q9T0Q9_BPFD/1-49 ['PDB; 1nh4 A; 1-49;']
COATB_BPIF1/22-73 ['PDB; 1ifk ; 1-50;']
for record in alignment:
    print(record)
ID: COATB_BPIKE/30-81
Name: COATB_BPIKE
Description: COATB_BPIKE/30-81
Database cross-references: PDB; 1ifl ; 1-52;
Number of features: 0
/accession=P03620.1
/start=30
/end=81
Per letter annotation for: secondary_structure
Seq('AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA')
ID: Q9T0Q8_BPIKE/1-52
Name: Q9T0Q8_BPIKE
Description: Q9T0Q8_BPIKE/1-52
Number of features: 0
/accession=Q9T0Q8.1
/start=1
/end=52
Seq('AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA')
ID: COATB_BPI22/32-83
Name: COATB_BPI22
Description: COATB_BPI22/32-83
Number of features: 0
/accession=P15416.1
/start=32
/end=83
Seq('DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA')
ID: COATB_BPM13/24-72
Name: COATB_BPM13
Description: COATB_BPM13/24-72
Database cross-references: PDB; 2cpb ; 1-49;, PDB; 2cps ; 1-49;
Number of features: 0
/accession=P69541.1
/start=24
/end=72
Per letter annotation for: secondary_structure
Seq('AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA')
ID: COATB_BPZJ2/1-49
Name: COATB_BPZJ2
Description: COATB_BPZJ2/1-49
Number of features: 0
/accession=P03618.1
/start=1
/end=49
Seq('AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA')
ID: Q9T0Q9_BPFD/1-49
Name: Q9T0Q9_BPFD
Description: Q9T0Q9_BPFD/1-49
Database cross-references: PDB; 1nh4 A; 1-49;
Number of features: 0
/accession=Q9T0Q9.1
/start=1
/end=49
Per letter annotation for: secondary_structure
Seq('AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA')
ID: COATB_BPIF1/22-73
Name: COATB_BPIF1
Description: COATB_BPIF1/22-73
Database cross-references: PDB; 1ifk ; 1-50;
Number of features: 0
/accession=P03619.2
/start=22
/end=73
Per letter annotation for: secondary_structure
Seq('FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA')
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord 
from Bio.Align import MultipleSeqAlignment 
align1 = MultipleSeqAlignment(
...     [
...          SeqRecord(Seq("ACTGCTAGCTAG"), id="Alpha"), 
...          SeqRecord(Seq("ACT-CTAGCTAG"), id="Beta"),
...          SeqRecord(Seq("ACTGCTAGDTAG"), id="Gamma"), 
...     ]
... )
>>> align2 = MultipleSeqAlignment(
...     [
...          SeqRecord(Seq("GTCAGC-AG"), id="Delta"), 
...          SeqRecord(Seq("GACAGCTAG"), id="Epsilon"),
...          SeqRecord(Seq("GTCAGCTAG"), id="Zeta"), 
...     ]
... )
>>> align3 = MultipleSeqAlignment(
...     [
...          SeqRecord(Seq("ACTAGTACAGCTG"), id="Eta"), 
...          SeqRecord(Seq("ACTAGTACAGCT-"), id="Theta"),
...          SeqRecord(Seq("-CTACTACAGGTG"), id="Iota"), 
...     ]
... )
my_alignments = [align1, align2, align3]
my_alignments
[<<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 12) at 7f4b6efdcf10>,
 <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 9) at 7f4b6efdc690>,
 <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 13) at 7f4b6efdc5d0>]
print(my_alignments)
[<<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 12) at 7f4b6efdcf10>, <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 9) at 7f4b6efdc690>, <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 13) at 7f4b6efdc5d0>]
from Bio import AlignIO
AlignIO.write(my_alignments, "my_example.phy", "phylip")
3
alignments = AlignIO.parse("my_example.phy", "phylip")
for alignment in alignments:
    print(alignment)
    print()
ID: Eta
Name: Eta
Description: Eta
Number of features: 0
Seq('ACTAGTACAGCTG')

ID: Theta
Name: Theta
Description: Theta
Number of features: 0
Seq('ACTAGTACAGCT-')

ID: Iota
Name: Iota
Description: Iota
Number of features: 0
Seq('-CTACTACAGGTG')

alignments = list(AlignIO.parse("my_example.phy", "phylip"))
last_align = alignments[-1]
print(last_align)
Alignment with 3 rows and 13 columns
ACTAGTACAGCTG Eta
ACTAGTACAGCT- Theta
-CTACTACAGGTG Iota
first_align = alignments[0]
print(first_align)
Alignment with 3 rows and 12 columns
ACTGCTAGCTAG Alpha
ACT-CTAGCTAG Beta
ACTGCTAGDTAG Gamma
from Bio import AlignIO
count = AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.aln", "clustal")
print ("Converted %i alignments" % count)
Converted 1 alignments
alignments = AlignIO.parse("PF05371_seed.sth", "stockholm")
count = AlignIO.write(alignments, "PF05371_seed.aln", "clustal" )
print("Converted %i alignments" % count)
Converted 1 alignments
AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.aln", "phylip" )
1
AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.phy", "phylip-relaxed")
1
alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
name_mapping = {}
for i,record in enumerate(alignment):
    name_mapping[i] = record.id
    record.id = "seq%i" %i
print(name_mapping)
{0: 'COATB_BPIKE/30-81', 1: 'Q9T0Q8_BPIKE/1-52', 2: 'COATB_BPI22/32-83', 3: 'COATB_BPM13/24-72', 4: 'COATB_BPZJ2/1-49', 5: 'Q9T0Q9_BPFD/1-49', 6: 'COATB_BPIF1/22-73'}
AlignIO.write([alignment], "PF0531_seed.phy", "phylip")
1
alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
print("Number of rows: %i" % len(alignment))
Number of rows: 7
for record in alignment:
    print("%s - %s" % (record.seq, record.id))
AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73
print(alignment[3:7])
Alignment with 4 rows and 52 columns
AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73
print(alignment[2, 6])
T
print(alignment[2].seq[6])
T
print(alignment[:,6])
TTT---T
print(alignment[3:6, :6])
Alignment with 3 rows and 6 columns
AEGDDP COATB_BPM13/24-72
AEGDDP COATB_BPZJ2/1-49
AEGDDP Q9T0Q9_BPFD/1-49
print(alignment[:, :6])
Alignment with 7 rows and 6 columns
AEPNAA COATB_BPIKE/30-81
AEPNAA Q9T0Q8_BPIKE/1-52
DGTSTA COATB_BPI22/32-83
AEGDDP COATB_BPM13/24-72
AEGDDP COATB_BPZJ2/1-49
AEGDDP Q9T0Q9_BPFD/1-49
FAADDA COATB_BPIF1/22-73
print(alignment[:, 6:9])
Alignment with 7 rows and 3 columns
TNY COATB_BPIKE/30-81
TNY Q9T0Q8_BPIKE/1-52
TSY COATB_BPI22/32-83
--- COATB_BPM13/24-72
--- COATB_BPZJ2/1-49
--- Q9T0Q9_BPFD/1-49
TSQ COATB_BPIF1/22-73
print(alignment[:, 9:])
Alignment with 7 rows and 43 columns
ATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
ATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
ATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
AKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73
edited = alignment[:, :6] + alignment[:, 9:]
print(edited)
Alignment with 7 rows and 49 columns
AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
DGTSTAATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
AEGDDPAKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
FAADDAAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73
edited.sort()
print(edited)
Alignment with 7 rows and 49 columns
DGTSTAATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
FAADDAAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73
AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
AEGDDPAKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord 
from Bio.Align import MultipleSeqAlignment
​
alignment = MultipleSeqAlignment(
    [
         SeqRecord(Seq("ACTCCTA"), id="seq1"),
         SeqRecord(Seq("AAT-CTA"), id="seq2"),
         SeqRecord(Seq("CCTACT-"), id="seq3"),
         SeqRecord(Seq("TCTCCTC"), id="seq4"),
    ]
)
print(alignment)
Alignment with 4 rows and 7 columns
ACTCCTA seq1
AAT-CTA seq2
CCTACT- seq3
TCTCCTC seq4
substitutions = alignment.substitutions
print(substitutions)
    A    C    T
A 2.0  4.5  1.0
C 4.5 10.0  0.5
T 1.0  0.5 12.0

m = substitutions.select("ATCG")
print(m)
    A    T    C   G
A 2.0  1.0  4.5 0.0
T 1.0 12.0  0.5 0.0
C 4.5  0.5 10.0 0.0
G 0.0  0.0  0.0 0.0

import Bio.Align.Applications
dir(Bio.Align.Applications)
['ClustalOmegaCommandline',
 'ClustalwCommandline',
 'DialignCommandline',
 'MSAProbsCommandline',
 'MafftCommandline',
 'MuscleCommandline',
 'PrankCommandline',
 'ProbconsCommandline',
 'TCoffeeCommandline',
 '_ClustalOmega',
 '_Clustalw',
 '_Dialign',
 '_MSAProbs',
 '_Mafft',
 '_Muscle',
 '_Prank',
 '_Probcons',
 '_TCoffee',
 '__all__',
 '__builtins__',
 '__cached__',
 '__doc__',
 '__file__',
 '__loader__',
 '__name__',
 '__package__',
 '__path__',
 '__spec__']
from Bio.Align.Applications import ClustalwCommandline
#https://raw.githubusercontent.com/biopython/biopython/refs/heads/master/Doc/examples/opuntia.aln
from Bio import AlignIO
align = AlignIO.read("opuntia.aln.txt", "clustal")
print(align)
Alignment with 7 rows and 906 columns
TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF191
TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF191
TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF191
TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF191
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF191
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF191
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273291|gb|AF191665.1|AF191
from Bio import Phylo
#https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/opuntia.dnd
tree = Phylo.read("opuntia.dnd.txt", "newick")
Phylo.draw_ascii(tree)
                             _______________ gi|6273291|gb|AF191665.1|AF191665
  __________________________|
 |                          |   ______ gi|6273290|gb|AF191664.1|AF191664
 |                          |__|
 |                             |_____ gi|6273289|gb|AF191663.1|AF191663
 |
_|_________________ gi|6273287|gb|AF191661.1|AF191661
 |
 |__________ gi|6273286|gb|AF191660.1|AF191660
 |
 |    __ gi|6273285|gb|AF191659.1|AF191659
 |___|
     | gi|6273284|gb|AF191658.1|AF191658

```
#Blast 

```


        from Bio.Blast import NCBIWWW
NCBIWWW.email = "dr.josh.vandenbrink@gmail.com"
result_handle = NCBIWWW.qblast("blastn", "nt", "8332116")
#https://github.com/biopython.biopython/blob/master/Doc/examples/m_cold.fasta
from Bio import SeqIO
record = SeqIO.read("m_cold.fasta.txt", format = "fasta")
print(record)
ID: gi|8332116|gb|BE037100.1|BE037100
Name: gi|8332116|gb|BE037100.1|BE037100
Description: gi|8332116|gb|BE037100.1|BE037100 MP14H09 MP Mesembryanthemum crystallinum cDNA 5' similar to cold acclimation protein, mRNA sequence
Number of features: 0
Seq('CACTAGTACTCGAGCGTNCTGCACCAATTCGGCACGAGCAAGTGACTACGTTNT...TTC')
result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)
with open("my_blast.xml", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()
from Bio.Blast import NCBIXML
result_handle = open("my_blast.xml")
blast_record = NCBIXML.read(result_handle)
E_VALUE_THRESH = 0.04
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print("****ALIGHTMENT****")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e_value:", hsp.expect)
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")
****ALIGHTMENT****
sequence: gi|1219041180|ref|XM_021875076.1| PREDICTED: Chenopodium quinoa cold-regulated 413 plasma membrane protein 2-like (LOC110697660), mRNA
length: 1173
e_value: 5.25852e-117
ACAGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTC...
|| ||||||||| |||| | |||| ||  |||| |||| | |||| ||| | |||| ||| ||| ||||| | ||...
ACCGAAAATGGGCAGAGGAGTGAATTATATGGCAATGACACCTGAGCAACTAGCCGCGGCCAATTTGATCAACTC...
****ALIGHTMENT****
sequence: gi|2514617377|ref|XM_021992092.2| PREDICTED: Spinacia oleracea cold-regulated 413 plasma membrane protein 2-like (LOC110787470), mRNA
length: 752
e_value: 1.41106e-111
AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
|||||||| |||  |||| | || ||||| |||||||| || ||||| |||| ||| ||| ||||||||||||||...
AAAATGGGTAGACGAATGGATTATTTGGCGATGAAAACCGAGCAATTAGCCGCGGCCAATTTGATCGATTCCGAT...
****ALIGHTMENT****
sequence: gi|2518612504|ref|XM_010682658.3| PREDICTED: Beta vulgaris subsp. vulgaris cold-regulated 413 plasma membrane protein 2 (LOC104895996), mRNA
length: 621
e_value: 3.78639e-106
TTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATATCAATGAGCTTAAAATGGCAACA...
||||||||||||||||| ||| ||||  |||||||| |||| ||||  ||||| ||||| ||||| || ||    ...
TTGGCCATGAAAACTGAGCAAATGGCGTTGGCTAATTTGATAGATTATGATATGAATGAACTTAAGATCGCTTTG...
****ALIGHTMENT****
sequence: gi|2031543140|ref|XM_041168865.1| PREDICTED: Juglans microcarpa x Juglans regia cold-regulated 413 plasma membrane protein 2-like (LOC121265293), mRNA
length: 1020
e_value: 1.32158e-105
AATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATAT...
||||||||| |||  | |  | |||||||||||||||||||    ||||  |||  || ||||||| || |||| ...
AATGGGGAG-GAA--GGATAATTTGGCCATGAAAACTGATCC---GGCCACGGCGGATTTGATCGACTCTGATAA...
****ALIGHTMENT****
sequence: gi|2618480339|ref|XM_048479995.2| PREDICTED: Ziziphus jujuba cold-regulated 413 plasma membrane protein 2 (LOC107424728), mRNA
length: 1028
e_value: 4.61277e-105
AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
|||||||||||    ||| |||  ||||| |||| |||||||| |   |||  |||| |  ||||  |||| |||...
AAAATGGGGAGG---ATGGAGTTTTTGGCTATGAGAACTGATCCA---GCCACGGCTGACTTGATAAATTCTGAT...
****ALIGHTMENT****
sequence: gi|2082357255|ref|XM_043119049.1| PREDICTED: Carya illinoinensis cold-regulated 413 plasma membrane protein 2-like (LOC122306609), transcript variant X2, mRNA
length: 1036
e_value: 5.6195e-104
ATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATATC...
|||||||| |||  | | || |||||||||||||||||||    ||||  |||  || ||||||| || ||||||...
ATGGGGAG-GAA--GGATTATTTGGCCATGAAAACTGATCC---GGCCACGGCGGATTTGATCGACTCTGATATC...
****ALIGHTMENT****
sequence: gi|2082357253|ref|XM_043119041.1| PREDICTED: Carya illinoinensis cold-regulated 413 plasma membrane protein 2-like (LOC122306609), transcript variant X1, mRNA
length: 1020
e_value: 5.6195e-104
ATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATATC...
|||||||| |||  | | || |||||||||||||||||||    ||||  |||  || ||||||| || ||||||...
ATGGGGAG-GAA--GGATTATTTGGCCATGAAAACTGATCC---GGCCACGGCGGATTTGATCGACTCTGATATC...
****ALIGHTMENT****
sequence: gi|1882610310|ref|XM_035691634.1| PREDICTED: Juglans regia cold-regulated 413 plasma membrane protein 2 (LOC108995251), transcript variant X2, mRNA
length: 909
e_value: 6.84595e-103
AATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATT---------GGCCGTGGCTAATATGATCGA...
||||||||| |||  | | || |||||||||||||||||||             ||||  |||  || |||||||...
AATGGGGAG-GAA--GGATTATTTGGCCATGAAAACTGATCCGGCCACGGCCACGGCCACGGCGGATTTGATCGA...
****ALIGHTMENT****
sequence: gi|1882610309|ref|XM_018970776.2| PREDICTED: Juglans regia cold-regulated 413 plasma membrane protein 2 (LOC108995251), transcript variant X1, mRNA
length: 1025
e_value: 6.84595e-103
AATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATT---------GGCCGTGGCTAATATGATCGA...
||||||||| |||  | | || |||||||||||||||||||             ||||  |||  || |||||||...
AATGGGGAG-GAA--GGATTATTTGGCCATGAAAACTGATCCGGCCACGGCCACGGCCACGGCGGATTTGATCGA...
****ALIGHTMENT****
sequence: gi|1350315638|ref|XM_006425719.2| PREDICTED: Citrus clementina cold-regulated 413 plasma membrane protein 2 (LOC18037141), transcript variant X3, mRNA
length: 893
e_value: 5.26316e-98
AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
|||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
****ALIGHTMENT****
sequence: gi|2395983798|ref|XM_006466623.4| PREDICTED: Citrus sinensis cold-regulated 413 plasma membrane protein 2 (LOC102620025), transcript variant X3, mRNA
length: 1052
e_value: 5.26316e-98
AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
|||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
****ALIGHTMENT****
sequence: gi|2395983796|ref|XM_025094967.2| PREDICTED: Citrus sinensis cold-regulated 413 plasma membrane protein 2 (LOC102620025), transcript variant X1, mRNA
length: 980
e_value: 5.26316e-98
AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
|||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
****ALIGHTMENT****
sequence: gi|1204884098|ref|XM_021445554.1| PREDICTED: Herrania umbratica cold-regulated 413 plasma membrane protein 2-like (LOC110429488), mRNA
length: 905
e_value: 5.26316e-98
AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
|||||||||||   ||| | || ||||| |||||||| ||||   | || |  |   || |||||  ||| ||||...
AAATGGGGAGA---ATGGACTATTTGGCTATGAAAACAGATCCTGTAGCAGAAG---ATTTGATCAGTTCTGATA...
****ALIGHTMENT****
sequence: gi|2395983800|ref|XM_006466626.4| PREDICTED: Citrus sinensis cold-regulated 413 plasma membrane protein 2 (LOC102620025), transcript variant X5, mRNA
length: 913
e_value: 5.26316e-98
AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
|||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
****ALIGHTMENT****
sequence: gi|1350315641|ref|XM_024180293.1| PREDICTED: Citrus clementina cold-regulated 413 plasma membrane protein 2 (LOC18037141), transcript variant X4, mRNA
length: 868
e_value: 5.26316e-98
AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
|||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
****ALIGHTMENT****
sequence: gi|1350315636|ref|XM_006425716.2| PREDICTED: Citrus clementina cold-regulated 413 plasma membrane protein 2 (LOC18037141), transcript variant X2, mRNA
length: 881
e_value: 5.26316e-98
AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
|||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
****ALIGHTMENT****
sequence: gi|1350315634|ref|XM_006425717.2| PREDICTED: Citrus clementina cold-regulated 413 plasma membrane protein 2 (LOC18037141), transcript variant X1, mRNA
length: 952
e_value: 5.26316e-98
AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
|||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
****ALIGHTMENT****
sequence: gi|2395983797|ref|XM_006466624.4| PREDICTED: Citrus sinensis cold-regulated 413 plasma membrane protein 2 (LOC102620025), transcript variant X2, mRNA
length: 968
e_value: 5.26316e-98
AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
|||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
****ALIGHTMENT****
sequence: gi|2395983799|ref|XM_006466625.3| PREDICTED: Citrus sinensis cold-regulated 413 plasma membrane protein 2 (LOC102620025), transcript variant X4, mRNA
length: 978
e_value: 5.26316e-98
AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
|||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
****ALIGHTMENT****
sequence: gi|1227938481|ref|XM_022049453.1| PREDICTED: Carica papaya cold-regulated 413 plasma membrane protein 2-like (LOC110820077), mRNA
length: 1009
e_value: 2.23795e-96
AGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCG...
|||||||||||||    ||| | || ||||| ||||| ||||||||   ||||   ||| || | |||  ||| |...
AGAAAATGGGGAGG---ATGGAATATTTGGCTATGAAGACTGATCA---GGCCACTGCTGATCTCATCACTTCTG...
****ALIGHTMENT****
sequence: gi|1063463253|ref|XM_007047033.2| PREDICTED: Theobroma cacao cold-regulated 413 plasma membrane protein 2 (LOC18611025), transcript variant X2, mRNA
length: 1071
e_value: 9.51602e-95
TGTGAACAGA-AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGAT...
|| |||| || |||||||||||   ||| | || ||||| |||||||| ||||   | || |  |   || ||||...
TGAGAACTGAGAAATGGGGAGA---ATGGACTATTTGGCTATGAAAACAGATCCTGTAGCAGAAG---ATTTGAT...
****ALIGHTMENT****
sequence: gi|1063463252|ref|XM_007047032.2| PREDICTED: Theobroma cacao cold-regulated 413 plasma membrane protein 2 (LOC18611025), transcript variant X1, mRNA
length: 1065
e_value: 9.51602e-95
TGTGAACAGA-AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGAT...
|| |||| || |||||||||||   ||| | || ||||| |||||||| ||||   | || |  |   || ||||...
TGAGAACTGAGAAATGGGGAGA---ATGGACTATTTGGCTATGAAAACAGATCCTGTAGCAGAAG---ATTTGAT...
****ALIGHTMENT****
sequence: gi|1269881403|ref|XM_022895603.1| PREDICTED: Durio zibethinus cold-regulated 413 plasma membrane protein 2 (LOC111300020), transcript variant X1, mRNA
length: 1072
e_value: 3.32142e-94
AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
|||||||||||   ||| |||| ||||| |||||||||||||   | || |  |  ||| |||||  ||| ||||...
AAATGGGGAGA---ATGGAGTATTTGGCTATGAAAACTGATCCTGTAGCTGAAG--AAT-TGATCAGTTCTGATA...
****ALIGHTMENT****
sequence: gi|1269881405|ref|XM_022895604.1| PREDICTED: Durio zibethinus cold-regulated 413 plasma membrane protein 2 (LOC111300020), transcript variant X2, mRNA
length: 1091
e_value: 3.32142e-94
AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
|||||||||||   ||| |||| ||||| |||||||||||||   | || |  |  ||| |||||  ||| ||||...
AAATGGGGAGA---ATGGAGTATTTGGCTATGAAAACTGATCCTGTAGCTGAAG--AAT-TGATCAGTTCTGATA...
****ALIGHTMENT****
sequence: gi|1269881407|ref|XM_022895605.1| PREDICTED: Durio zibethinus cold-regulated 413 plasma membrane protein 2 (LOC111300020), transcript variant X3, mRNA
length: 1069
e_value: 3.32142e-94
AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
|||||||||||   ||| |||| ||||| |||||||||||||   | || |  |  ||| |||||  ||| ||||...
AAATGGGGAGA---ATGGAGTATTTGGCTATGAAAACTGATCCTGTAGCTGAAG--AAT-TGATCAGTTCTGATA...
****ALIGHTMENT****
sequence: gi|2082386143|ref|XM_043113301.1| PREDICTED: Carya illinoinensis cold-regulated 413 plasma membrane protein 2-like (LOC122301958), transcript variant X1, mRNA
length: 844
e_value: 1.15929e-93
ATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATATCAATGAGCTTAAA...
||||| |||||||| |||||||||||||    |||  |||   || ||||||  || ||||||||||| || || ...
ATGAATTACTTGGCTATGAAAACTGATCC---GGCAATGGAGGATTTGATCGGCTCTGATATCAATGACCTCAAG...
****ALIGHTMENT****
sequence: gi|2082386146|ref|XM_043113302.1| PREDICTED: Carya illinoinensis cold-regulated 413 plasma membrane protein 2-like (LOC122301958), transcript variant X2, mRNA
length: 824
e_value: 1.15929e-93
ATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATATCAATGAGCTTAAA...
||||| |||||||| |||||||||||||    |||  |||   || ||||||  || ||||||||||| || || ...
ATGAATTACTTGGCTATGAAAACTGATCC---GGCAATGGAGGATTTGATCGGCTCTGATATCAATGACCTCAAG...
****ALIGHTMENT****
sequence: gi|1954740698|ref|XM_038867092.1| PREDICTED: Tripterygium wilfordii cold-regulated 413 plasma membrane protein 2 (LOC120014952), mRNA
length: 999
e_value: 4.04632e-93
GAACAGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGAT...
||| ||||||||||||||   | | | || ||||| ||||| |||||||    ||  ||||   || |||||   ...
GAAAAGAAAATGGGGAGA---ACGGATTATTTGGCGATGAAGACTGATCC---GGTTGTGGACGATTTGATCAGC...
****ALIGHTMENT****
sequence: gi|1882636119|ref|XM_018974650.2| PREDICTED: Juglans regia cold-regulated 413 plasma membrane protein 2-like (LOC108998174), mRNA
length: 1015
e_value: 4.92942e-92
AATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATAT...
|||||||||    ||||| || ||||| |||||||||||||    |||  ||| | || ||||||  || |||||...
AATGGGGAGG---ATGAATTATTTGGCTATGAAAACTGATCC---GGCAATGGATGATTTGATCGGCTCTGATAT...
****ALIGHTMENT****
sequence: gi|2526866810|ref|XM_057645500.1| PREDICTED: Actinidia eriantha cold-regulated 413 plasma membrane protein 2-like (LOC130785340), mRNA
length: 1152
e_value: 4.92942e-92
AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
||||||||||||   ||| | || ||||| ||||| || |||| |  |||  || | ||| ||||| |||| || ...
AAAATGGGGAGA---ATGGATTATTTGGCGATGAAGACCGATCCAGCGGC--TGCCGAAT-TGATCAATTCGGAC...
****ALIGHTMENT****
sequence: gi|1187397285|gb|KX009413.1| Santalum album COR413-PM2 mRNA, complete cds
length: 837
e_value: 1.72054e-91
AATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATAT...
|||||||||    ||| | | ||||||||||||||| ||||    |||||  |   || ||||| ||||||| ||...
AATGGGGAGG---ATGGATTTCTTGGCCATGAAAACAGATCCCGCGGCCGCCG---ATTTGATCAATTCCGACAT...
****ALIGHTMENT****
sequence: gi|2550782781|ref|XM_058372567.1| PREDICTED: Rhododendron vialii cold-regulated 413 plasma membrane protein 2 (LOC131336659), mRNA
length: 1110
e_value: 2.09604e-90
GCCGTGGCTAATATGATCGATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGT...
||||  ||| | |||||||| || |||||||| ||||| || || ||  | | | | || || |   | || |  ...
GCCGATGCTGAAATGATCGACTCGGATATCAACGAGCTGAAGATCGCGGCCAAGCGACTGATTAGCCACGCCACC...
****ALIGHTMENT****
sequence: gi|2806124758|ref|XM_068481225.1| PREDICTED: Pyrus communis cold-regulated 413 plasma membrane protein 2-like (LOC137741519), mRNA
length: 850
e_value: 7.31591e-90
TGATCGATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGTATGCTCGGTCATT...
|||| ||||| |||||||| ||||| || || ||| | | ||| ||||||| |||||| |  | ||| |||  ||...
TGATAGATTCAGATATCAAAGAGCTCAAGATTGCAGCCAAGAGACTCATCAGTGATGCCACCAAGCTTGGTGGTT...
****ALIGHTMENT****
sequence: gi|2532162279|ref|XM_058104265.1| PREDICTED: Malania oleifera cold-regulated 413 plasma membrane protein 2-like (LOC131152402), mRNA
length: 2364
e_value: 8.9126e-89
GAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGA...
||||||||||||      | |||||||||||||||||||||||| |  ||| |  |   || |||||  ||| ||...
GAAAATGGGGAGGTC---GGAGTACTTGGCCATGAAAACTGATCCAGCGGCTGCCG---ATTTGATCAGTTCGGA...
****ALIGHTMENT****
sequence: gi|2250518185|ref|XM_009343631.3| PREDICTED: Pyrus x bretschneideri cold-regulated 413 plasma membrane protein 2 (LOC103933927), mRNA
length: 787
e_value: 8.9126e-89
TGATCGATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGTATGCTCGGTCATT...
|||| ||||| |||||||| ||||| || || ||| | | ||| ||||||| |||||| |  | ||| |||  ||...
TGATAGATTCAGATATCAAAGAGCTCAAGATTGCAGCCAAGAGACTCATCAGTGATGCCACCAAGCTTGGTGGTT...
****ALIGHTMENT****
sequence: gi|1350280614|ref|XM_024170292.1| PREDICTED: Morus notabilis cold-regulated 413 plasma membrane protein 2 (LOC21394987), mRNA
length: 1020
e_value: 3.1108e-88
AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
|||||||||| ||       || |||||||||||||| || | |   |||  |||| || ||||  |||| ||||...
AAATGGGGAGGGAT------TATTTGGCCATGAAAACGGACCCA---GCCACGGCTGATTTGATAAATTCTGATA...
****ALIGHTMENT****
sequence: gi|743838297|ref|XM_011027373.1| PREDICTED: Populus euphratica cold-regulated 413 plasma membrane protein 2 (LOC105126500), transcript variant X2, mRNA
length: 1132
e_value: 3.1108e-88
AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
|||||||||||    ||| |||  |||   ||||| |||||| |    | |   |||||| | || |||||||||...
AAAATGGGGAGG---ATGGAGTTTTTGAAGATGAAGACTGATGATGAAGTCAGCGCTAATTTAATTGATTCCGAT...
****ALIGHTMENT****
sequence: gi|743838293|ref|XM_011027372.1| PREDICTED: Populus euphratica cold-regulated 413 plasma membrane protein 2 (LOC105126500), transcript variant X1, mRNA
length: 980
e_value: 3.1108e-88
AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
|||||||||||    ||| |||  |||   ||||| |||||| |    | |   |||||| | || |||||||||...
AAAATGGGGAGG---ATGGAGTTTTTGAAGATGAAGACTGATGATGAAGTCAGCGCTAATTTAATTGATTCCGAT...
****ALIGHTMENT****
sequence: gi|1768569081|ref|XM_031406607.1| PREDICTED: Pistacia vera cold-regulated 413 plasma membrane protein 2-like (LOC116120644), mRNA
length: 982
e_value: 3.78974e-87
AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATT-GGCCGTGGCTAATATGATCGATTCCGA...
|||||||||||    ||| | ||  |||  ||||||||||| ||||  ||     ||| |  ||||  | || ||...
AAAATGGGGAGG---ATGGATTATCTGGGAATGAAAACTGA-CAATCAGGTTACTGCTGAGGTGATTAACTCTGA...
****ALIGHTMENT****
sequence: gi|2396494064|ref|XM_024605027.2| PREDICTED: Populus trichocarpa cold-regulated 413 plasma membrane protein 2 (LOC18101203), transcript variant X2, mRNA
length: 1178
e_value: 1.32275e-86
AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
|||||||||||    ||| |||  |||   ||||| |||||| |    | |   |||||| | || || ||||||...
AAAATGGGGAGG---ATGGAGTTTTTGAAGATGAAGACTGATGATGAAGTCAGCGCTAATTTAATTGAGTCCGAT...
****ALIGHTMENT****
sequence: gi|2396494060|ref|XM_052454347.1| PREDICTED: Populus trichocarpa cold-regulated 413 plasma membrane protein 2 (LOC18101203), transcript variant X1, mRNA
length: 1018
e_value: 1.32275e-86
AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
|||||||||||    ||| |||  |||   ||||| |||||| |    | |   |||||| | || || ||||||...
AAAATGGGGAGG---ATGGAGTTTTTGAAGATGAAGACTGATGATGAAGTCAGCGCTAATTTAATTGAGTCCGAT...
****ALIGHTMENT****
sequence: gi|1585724761|ref|XM_028202722.1| PREDICTED: Camellia sinensis cold-regulated 413 plasma membrane protein 2-like (LOC114262355), mRNA
length: 910
e_value: 4.61684e-86
AGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCG...
|||||||||||||  ||||| |||| ||||| ||||| || |||||    |||    |  |   |||  ||||||...
AGAAAATGGGGAGGAAAATGGAGTATTTGGCAATGAAGACCGATCATCCAGCCCCAACCCAATCGATGAATTCCG...
****ALIGHTMENT****
sequence: gi|2537663858|ref|XM_021815584.2| PREDICTED: Hevea brasiliensis cold-regulated 413 plasma membrane protein 2 (LOC110658100), mRNA
length: 945
e_value: 1.61144e-85
AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
||||||||||    ||| ||||||||   ||||  ||||||||| |  |   |||  || ||||| | || ||| ...
AAATGGGGAGG---ATGGAGTACTTGAAAATGAGTACTGATCAAGTACC---GGCCGATTTGATCAAGTCTGATC...
****ALIGHTMENT****
sequence: gi|1860377401|ref|XM_035077206.1| PREDICTED: Populus alba cold-regulated 413 plasma membrane protein 2-like (LOC118063227), transcript variant X2, mRNA
length: 916
e_value: 1.61144e-85
AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
|||||||||||    ||| |||  |||   ||||| |||||| |    | |   | |||| | || || ||||||...
AAAATGGGGAGG---ATGGAGTTTTTGAAGATGAAGACTGATGATGAAGTCAGCGGTAATTTAATTGAGTCCGAT...
****ALIGHTMENT****
sequence: gi|2645357626|ref|XM_062094449.1| PREDICTED: Populus nigra cold-regulated 413 plasma membrane protein 2-like (LOC133673573), mRNA
length: 1175
e_value: 1.61144e-85
AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
|||||||||||    ||| |||  |||   ||||| |||||| |    | |   |||||| | || || ||||||...
AAAATGGGGAGG---ATGGAGTTTTTGAAGATGAAGACTGATGATGAAGTCAGCGCTAATTTAATTGAGTCCGAT...
****ALIGHTMENT****
sequence: gi|1860377399|ref|XM_035077205.1| PREDICTED: Populus alba cold-regulated 413 plasma membrane protein 2-like (LOC118063227), transcript variant X1, mRNA
length: 1109
e_value: 1.61144e-85
AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
|||||||||||    ||| |||  |||   ||||| |||||| |    | |   | |||| | || || ||||||...
AAAATGGGGAGG---ATGGAGTTTTTGAAGATGAAGACTGATGATGAAGTCAGCGGTAATTTAATTGAGTCCGAT...
****ALIGHTMENT****
sequence: gi|1162571918|ref|XM_007202530.2| PREDICTED: Prunus persica cold-regulated 413 plasma membrane protein 2 (LOC18770198), transcript variant X1, mRNA
length: 811
e_value: 1.96313e-84
TGATCGATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGTATGCTCGGTCATT...
||||  |||| || |||||||| || || || ||| | | ||  |||||||||||||| | || ||| |||   |...
TGATAAATTCAGACATCAATGATCTCAAGATTGCAGCCAAGAAACTCATCAATGATGCCACTAAGCTTGGTGGGT...
****ALIGHTMENT****
sequence: gi|1162571919|ref|XM_020568695.1| PREDICTED: Prunus persica cold-regulated 413 plasma membrane protein 2 (LOC18770198), transcript variant X2, mRNA
length: 929
e_value: 1.96313e-84
TGATCGATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGTATGCTCGGTCATT...
||||  |||| || |||||||| || || || ||| | | ||  |||||||||||||| | || ||| |||   |...
TGATAAATTCAGACATCAATGATCTCAAGATTGCAGCCAAGAAACTCATCAATGATGCCACTAAGCTTGGTGGGT...
****ALIGHTMENT****
sequence: gi|2583747300|ref|XM_059787294.1| PREDICTED: Cornus florida cold-regulated 413 plasma membrane protein 2-like (LOC132285128), mRNA
length: 1126
e_value: 1.96313e-84
AGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCG...
|||||||||||||| |   | |||| ||||| |||||||||||||    ||||   ||  |  ||||| ||||||...
AGAAAATGGGGAGAAA---GGAGTATTTGGCTATGAAAACTGATCC---GGCCACAGCCGAATTGATCAATTCCG...
****ALIGHTMENT****
sequence: gi|1229761331|ref|XM_022277554.1| PREDICTED: Momordica charantia cold-regulated 413 plasma membrane protein 2-like (LOC111005887), mRNA
length: 850
e_value: 6.852e-84
ATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGTATGCTCGGTCATTACGGGT...
|||| |||||||| ||||||||||| ||| | | ||||||  |  |  |||| |  | |||||||    | ||  ...
ATTCTGATATCAACGAGCTTAAAATTGCAGCCACGAGGCTTCTTGAACATGCCACCAAGCTCGGTGGAAAGGGCC...


```


# Challenge 1

```



Code

Python 3
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIWWW
NCBIWWW.email = "gqu002@email.latech.edu"
result_handle = NCBIWWW.qblast("blastn", "nt", "37455")
from Bio import SeqIO
record = SeqIO.read("sequence.fasta", format = "fasta")
print(record)
ID: NT_033778.4:21522420-21559977
Name: NT_033778.4:21522420-21559977
Description: NT_033778.4:21522420-21559977 Drosophila melanogaster chromosome 2R
Number of features: 0
Seq('TAAATTTGAAGATAAATTTCAAGCCGGGCTGCCCGTGTTTCAGTTTCCCAACGA...TTC')
result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)
with open("my_blast.xml", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()
from Bio.Blast import NCBIXML
result_handle = open("my_blast.xml")
blast_record = NCBIXML.read(result_handle)
E_VALUE_THRESH = 0.04
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print("****ALIGHTMENT****")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e_value:", hsp.expect)
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")
****ALIGHTMENT****
sequence: gi|16874826|gb|AC099307.1| Drosophila melanogaster, chromosome 2R, region 57C-57D, BAC clone BACR13H23, complete sequence
length: 179139
e_value: 0.0
TAAATTTGAAGATAAATTTCAAGCCGGGCTGCCCGTGTTTCAGTTTCCCAACGAAAACATCAACAAATAACATCA...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
TAAATTTGAAGATAAATTTCAAGCCGGGCTGCCCGTGTTTCAGTTTCCCAACGAAAACATCAACAAATAACATCA...
****ALIGHTMENT****
sequence: gi|3687207|gb|AC004300.1|AC004300 Drosophila melanogaster DNA sequence (P1 DS08209 (D186)), complete sequence
length: 58663
e_value: 0.0
TAAATTTGAAGATAAATTTCAAGCCGGGCTGCCCGTGTTTCAGTTTCCCAACGAAAACATCAACAAATAACATCA...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
TAAATTTGAAGATAAATTTCAAGCCGGGCTGCCCGTGTTTCAGTTTCCCAACGAAAACATCAACAAATAACATCA...
****ALIGHTMENT****
sequence: gi|15451509|gb|AC009459.6| Drosophila melanogaster, chromosome 2R, region 57E-57F, BAC clone BACR34C17, complete sequence
length: 167447
e_value: 0.0
GAATTCCAAGGCACAGGGCTATGACCATGCCAAAGGGCCAAATAAACAGGTAATCCAATTAAGCGGCTTCAGGTG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GAATTCCAAGGCACAGGGCTATGACCATGCCAAAGGGCCAAATAAACAGGTAATCCAATTAAGCGGCTTCAGGTG...
****ALIGHTMENT****
sequence: gi|40365018|gb|AY461301.1| Drosophila melanogaster strain NC086 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7400
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40365092|gb|AY461338.1| Drosophila melanogaster strain NC129 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7403
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40364933|gb|AY461258.1| Drosophila melanogaster strain NC041 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7399
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40364817|gb|AY461200.1| Drosophila melanogaster strain K3756 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7418
e_value: 0.0
CAGCAATAGCAATGCTCGATAATTTGATCGAATGGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATC...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
CAGCAATAGCAATGCTCGATAATTTGATCGAATGGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATC...
****ALIGHTMENT****
sequence: gi|40364811|gb|AY461197.1| Drosophila melanogaster strain K3748 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7418
e_value: 0.0
CAGCAATAGCAATGCTCGATAATTTGATCGAATGGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATC...
|||||||||||||||||||||||||||||||||||| |||||||||||| |||||||||||||||||||||||||...
CAGCAATAGCAATGCTCGATAATTTGATCGAATGGGATAACCATTCGAGCTCATACAGATACACGGCTCCAAATC...
****ALIGHTMENT****
sequence: gi|40364819|gb|AY461201.1| Drosophila melanogaster strain K3758 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7422
e_value: 0.0
CAGCAATAGCAATGCTCGATAATTTGATCGAATGGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATC...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
CAGCAATAGCAATGCTCGATAATTTGATCGAATGGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATC...
****ALIGHTMENT****
sequence: gi|40365060|gb|AY461322.1| Drosophila melanogaster strain NC111 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7385
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40364823|gb|AY461203.1| Drosophila melanogaster strain K3683 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7416
e_value: 0.0
CAGCAATAGCAATGCTCGATAATTTGATCGAATGGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATC...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
CAGCAATAGCAATGCTCGATAATTTGATCGAATGGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATC...
****ALIGHTMENT****
sequence: gi|40365022|gb|AY461303.1| Drosophila melanogaster strain NC088 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7398
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
|||||||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGGTAACCATTCGAGCTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40364877|gb|AY461230.1| Drosophila melanogaster strain NC011 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7390
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
|||||||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGGTAACCATTCGAGCTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40364863|gb|AY461223.1| Drosophila melanogaster strain NC003 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7408
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
||| |||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGATAACCATTCGAGCTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40364982|gb|AY461283.1| Drosophila melanogaster strain NC066 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7397
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40364883|gb|AY461233.1| Drosophila melanogaster strain NC014 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7392
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
|||||||||||||||| |||||||||||||||||||||||||||||||| |||||||||||||||||||||||||...
GGGGTAACCATTCGAGCTCATACAGATACACGGCTCCAAATCTCAAACGATGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40365128|gb|AY461356.1| Drosophila melanogaster strain NC149 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7407
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
||| |||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGATAACCATTCGAGCTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40365096|gb|AY461340.1| Drosophila melanogaster strain NC131 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7404
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40364945|gb|AY461264.1| Drosophila melanogaster strain NC047 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7311
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
|||||||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGGTAACCATTCGAGCTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40364927|gb|AY461255.1| Drosophila melanogaster strain NC038 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7408
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40365030|gb|AY461307.1| Drosophila melanogaster strain NC094 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7407
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40364931|gb|AY461257.1| Drosophila melanogaster strain NC040 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7404
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40364998|gb|AY461291.1| Drosophila melanogaster strain NC075 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7405
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40364957|gb|AY461270.1| Drosophila melanogaster strain NC053 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7392
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40364825|gb|AY461204.1| Drosophila melanogaster strain K3684 epidermal growth factor receptor (Egfr) gene, partial cds >gi|40364827|gb|AY461205.1| Drosophila melanogaster strain K3685 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7414
e_value: 0.0
CAGCAATAGCAATGCTCGATAATTTGATCGAATGGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATC...
|||||||||||||||| |||||||||||||||||||||||||||||||| |||||||||||||||||||||||||...
CAGCAATAGCAATGCTTGATAATTTGATCGAATGGGGTAACCATTCGAGCTCATACAGATACACGGCTCCAAATC...
****ALIGHTMENT****
sequence: gi|40365056|gb|AY461320.1| Drosophila melanogaster strain NC108 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7403
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40364831|gb|AY461207.1| Drosophila melanogaster strain K3692 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7419
e_value: 0.0
CAGCAATAGCAATGCTCGATAATTTGATCGAATGGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATC...
||||||||||||||||||||||||||||||||||||||||||||||||| |||||||||||||||||||||||||...
CAGCAATAGCAATGCTCGATAATTTGATCGAATGGGGTAACCATTCGAGCTCATACAGATACACGGCTCCAAATC...
****ALIGHTMENT****
sequence: gi|40365004|gb|AY461294.1| Drosophila melanogaster strain NC079 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7406
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40364976|gb|AY461280.1| Drosophila melanogaster strain NC063 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7305
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40364935|gb|AY461259.1| Drosophila melanogaster strain NC042 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7356
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40365050|gb|AY461317.1| Drosophila melanogaster strain NC105 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7396
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40364994|gb|AY461289.1| Drosophila melanogaster strain NC073 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7354
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40365014|gb|AY461299.1| Drosophila melanogaster strain NC084 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7402
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40365014|gb|AY461299.1| Drosophila melanogaster strain NC084 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7402
e_value: 0.0
CGCCATGCCAAATTCGAAGAGATCCCGGATCCGGATCGGCTCCGCCGGAGGGCACACACTAAGCTAAGCCGATGT...
||||||||||||||||||||||||||||||||||||||||||||||||||||| |||||||||||||||||||||...
CGCCATGCCAAATTCGAAGAGATCCCGGATCCGGATCGGCTCCGCCGGAGGGCGCACACTAAGCTAAGCCGATGT...
****ALIGHTMENT****
sequence: gi|40364861|gb|AY461222.1| Drosophila melanogaster strain NC002 epidermal growth factor receptor (Egfr) gene, partial cds
length: 6710
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
||| |||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGATAACCATTCGAGCTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40364984|gb|AY461284.1| Drosophila melanogaster strain NC068 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7393
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
||| |||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGATAACCATTCGAGCTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40364984|gb|AY461284.1| Drosophila melanogaster strain NC068 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7393
e_value: 0.0
TTAGACGCCTCCGCCATGCCAAATTCGAAGAGATCCCGGATCCGGATCGGCTCCGCCGGAGGGCACACACTAAGC...
|||||||||||||||||||||||||||||||||||||||||||||||||||||| ||||||||| ||||||||||...
TTAGACGCCTCCGCCATGCCAAATTCGAAGAGATCCCGGATCCGGATCGGCTCCACCGGAGGGCGCACACTAAGC...
****ALIGHTMENT****
sequence: gi|40365112|gb|AY461348.1| Drosophila melanogaster strain NC139 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7405
e_value: 0.0
TAATCAGTACGTCAGATGATGATCACCGTCGTGATGATGGCCAGAGGGGCGGCGAAGGGCGAATGGAGTCGTGTG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
TAATCAGTACGTCAGATGATGATCACCGTCGTGATGATGGCCAGAGGGGCGGCGAAGGGCGAATGGAGTCGTGTG...
****ALIGHTMENT****
sequence: gi|40365112|gb|AY461348.1| Drosophila melanogaster strain NC139 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7405
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
||| |||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGATAACCATTCGAGCTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40364885|gb|AY461234.1| Drosophila melanogaster strain NC015 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7408
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40364885|gb|AY461234.1| Drosophila melanogaster strain NC015 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7408
e_value: 0.0
CTTCTCTCGGGCGTCAATTACTTCTTGTACTTGTAGTTAGTAGTAGTCGTCTCTTTTAGCGCCTAGTTTCCAACT...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
CTTCTCTCGGGCGTCAATTACTTCTTGTACTTGTAGTTAGTAGTAGTCGTCTCTTTTAGCGCCTAGTTTCCAACT...
****ALIGHTMENT****
sequence: gi|40365100|gb|AY461342.1| Drosophila melanogaster strain NC133 epidermal growth factor receptor (Egfr) gene, partial cds
length: 6718
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
||| |||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGATAACCATTCGAGCTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40365102|gb|AY461343.1| Drosophila melanogaster strain NC134 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7407
e_value: 0.0
TGATGATCACCGTCGTGATGATGGCCAGAGGGGCGGCGAAGGGCGAATGGAGTCGTGTGTGTTGATGGCTTGGCA...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
TGATGATCACCGTCGTGATGATGGCCAGAGGGGCGGCGAAGGGCGAATGGAGTCGTGTGTGTTGATGGCTTGGCA...
****ALIGHTMENT****
sequence: gi|40365102|gb|AY461343.1| Drosophila melanogaster strain NC134 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7407
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
||| |||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGATAACCATTCGAGCTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40365024|gb|AY461304.1| Drosophila melanogaster strain NC089 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7407
e_value: 0.0
TAATCAGTACGTCAGATGATGATCACCGTCGTGATGATGGCCAGAGGGGCGGCGAAGGGCGAATGGAGTCGTGTG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
TAATCAGTACGTCAGATGATGATCACCGTCGTGATGATGGCCAGAGGGGCGGCGAAGGGCGAATGGAGTCGTGTG...
****ALIGHTMENT****
sequence: gi|40365024|gb|AY461304.1| Drosophila melanogaster strain NC089 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7407
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40364815|gb|AY461199.1| Drosophila melanogaster strain K3751 epidermal growth factor receptor (Egfr) gene, partial cds
length: 6689
e_value: 0.0
CACCGTCGTGATGATGGCCAGAGGGGCGGCGAAGGGCGAATGGAGTCGTGTGTGTTGATGGCTTGGCATACAAAT...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
CACCGTCGTGATGATGGCCAGAGGGGCGGCGAAGGGCGAATGGAGTCGTGTGTGTTGATGGCTTGGCATACAAAT...
****ALIGHTMENT****
sequence: gi|40364965|gb|AY461274.1| Drosophila melanogaster strain NC057 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7402
e_value: 0.0
ATGCCCAGCCATGCGCTCGGTTCGCAAGAAACCCGATGCTGATGCTGAAGCGTTTGCCTGGCAGGGAATCCACCG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
ATGCCCAGCCATGCGCTCGGTTCGCAAGAAACCCGATGCTGATGCTGAAGCGTTTGCCTGGCAGGGAATCCACCG...
****ALIGHTMENT****
sequence: gi|40364965|gb|AY461274.1| Drosophila melanogaster strain NC057 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7402
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40365120|gb|AY461352.1| Drosophila melanogaster strain NC144 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7399
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
||| |||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGATAACCATTCGAGCTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40365120|gb|AY461352.1| Drosophila melanogaster strain NC144 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7399
e_value: 0.0
AAATTGGCTCCCACCACCGATGGGTCCGAAGCCATTGCGGAACCCGATGACTACCTGCAACCCAAGGCAGCACCT...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
AAATTGGCTCCCACCACCGATGGGTCCGAAGCCATTGCGGAACCCGATGACTACCTGCAACCCAAGGCAGCACCT...
****ALIGHTMENT****
sequence: gi|40364833|gb|AY461208.1| Drosophila melanogaster strain K3697 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7426
e_value: 0.0
AAATAAAACGACCAGAAACAAAGCTTAAGAATTTGTACGAAAATTTTTGCGCTTCGATTGCTGGAGCCGCTTTTG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
AAATAAAACGACCAGAAACAAAGCTTAAGAATTTGTACGAAAATTTTTGCGCTTCGATTGCTGGAGCCGCTTTTG...
****ALIGHTMENT****
sequence: gi|40364833|gb|AY461208.1| Drosophila melanogaster strain K3697 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7426
e_value: 0.0
CAGCAATAGCAATGCTCGATAATTTGATCGAATGGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATC...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
CAGCAATAGCAATGCTCGATAATTTGATCGAATGGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATC...
****ALIGHTMENT****
sequence: gi|40364867|gb|AY461225.1| Drosophila melanogaster strain NC005 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7398
e_value: 0.0
ATGATCACCGTCGTGATGATGGCCAGAGGGGCGGCGAAGGGCGAATGGAGTCGTGTGTGTTGATGGCTTGGCATA...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
ATGATCACCGTCGTGATGATGGCCAGAGGGGCGGCGAAGGGCGAATGGAGTCGTGTGTGTTGATGGCTTGGCATA...
****ALIGHTMENT****
sequence: gi|40364867|gb|AY461225.1| Drosophila melanogaster strain NC005 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7398
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
||| |||||||||||| |||||||| | |||||||||||||||||||||||||||||||||||||||||||||||...
GGGATAACCATTCGAGCTCATACAGNTNCACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40364867|gb|AY461225.1| Drosophila melanogaster strain NC005 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7398
e_value: 0.0
TGTAGAATCCATGTAATTACGAGAACAAAACACCGATAGTGCATTTCTTAGACGCCTCCGCCATGCCAAATTCGA...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
TGTAGAATCCATGTAATTACGAGAACAAAACACCGATAGTGCATTTCTTAGACGCCTCCGCCATGCCAAATTCGA...
****ALIGHTMENT****
sequence: gi|40365046|gb|AY461315.1| Drosophila melanogaster strain NC103 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7412
e_value: 0.0
ATAAATAAAACGACCAGAAACAAAGCTTAAGAATTTGTACGAAAATTTTTGCGCTTCGATTGCTGGAGCCGCTTT...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
ATAAATAAAACGACCAGAAACAAAGCTTAAGAATTTGTACGAAAATTTTTGCGCTTCGATTGCTGGAGCCGCTTT...
****ALIGHTMENT****
sequence: gi|40365046|gb|AY461315.1| Drosophila melanogaster strain NC103 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7412
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
|||||||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGGTAACCATTCGAGCTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40364949|gb|AY461266.1| Drosophila melanogaster strain NC049 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7398
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40364949|gb|AY461266.1| Drosophila melanogaster strain NC049 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7398
e_value: 0.0
TGAGAAGGATCTCATCCGAAAATTGGCTCCCACCACCGATGGGTCCGAAGCCATTGCGGAACCCGATGACTACCT...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
TGAGAAGGATCTCATCCGAAAATTGGCTCCCACCACCGATGGGTCCGAAGCCATTGCGGAACCCGATGACTACCT...
****ALIGHTMENT****
sequence: gi|40364869|gb|AY461226.1| Drosophila melanogaster strain NC006 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7382
e_value: 0.0
ATTATGGGATATTTGTTCGAAGGAAATCGATAAATTAGGGAAGGAGGAGCGGCCATTTTGTGGATAATTAGTTGG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
ATTATGGGATATTTGTTCGAAGGAAATCGATAAATTAGGGAAGGAGGAGCGGCCATTTTGTGGATAATTAGTTGG...
****ALIGHTMENT****
sequence: gi|40364869|gb|AY461226.1| Drosophila melanogaster strain NC006 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7382
e_value: 0.0
ATGATGATCACCGTCGTGATGATGGCCAGAGGGGCGGCGAAGGGCGAATGGAGTCGTGTGTGTTGATGGCTTGGC...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
ATGATGATCACCGTCGTGATGATGGCCAGAGGGGCGGCGAAGGGCGAATGGAGTCGTGTGTGTTGATGGCTTGGC...
****ALIGHTMENT****
sequence: gi|40364869|gb|AY461226.1| Drosophila melanogaster strain NC006 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7382
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
||| |||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGATAACCATTCGAGCTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40364873|gb|AY461228.1| Drosophila melanogaster strain NC008 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7399
e_value: 0.0
ATTATGGGATATTTGTTCGAAGGAAATCGATAAATTAGGGAAGGAGGAGCGGCCATTTTGTGGATAATTAGTTGG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
ATTATGGGATATTTGTTCGAAGGAAATCGATAAATTAGGGAAGGAGGAGCGGCCATTTTGTGGATAATTAGTTGG...
****ALIGHTMENT****
sequence: gi|40364873|gb|AY461228.1| Drosophila melanogaster strain NC008 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7399
e_value: 0.0
ATGATCACCGTCGTGATGATGGCCAGAGGGGCGGCGAAGGGCGAATGGAGTCGTGTGTGTTGATGGCTTGGCATA...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
ATGATCACCGTCGTGATGATGGCCAGAGGGGCGGCGAAGGGCGAATGGAGTCGTGTGTGTTGATGGCTTGGCATA...
****ALIGHTMENT****
sequence: gi|40364873|gb|AY461228.1| Drosophila melanogaster strain NC008 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7399
e_value: 0.0
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
****ALIGHTMENT****
sequence: gi|40365040|gb|AY461312.1| Drosophila melanogaster strain NC100 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7394
e_value: 0.0
TAATCAGTACGTCAGATGATGATCACCGTCGTGATGATGGCCAGAGGGGCGGCGAAGGGCGAATGGAGTCGTGTG...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
TAATCAGTACGTCAGATGATGATCACCGTCGTGATGATGGCCAGAGGGGCGGCGAAGGGCGAATGGAGTCGTGTG...
****ALIGHTMENT****
sequence: gi|40365040|gb|AY461312.1| Drosophila melanogaster strain NC100 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7394
e_value: 0.0
CTAGCAACAAGAATTCGAGTACCGGAGACGATGAGACGGATTCGAGTGCCCGGGAAGTGGGCGTGGGTAATCTGC...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
CTAGCAACAAGAATTCGAGTACCGGAGACGATGAGACGGATTCGAGTGCCCGGGAAGTGGGCGTGGGTAATCTGC...
****ALIGHTMENT****
sequence: gi|40365040|gb|AY461312.1| Drosophila melanogaster strain NC100 epidermal growth factor receptor (Egfr) gene, partial cds
length: 7394
e_value: 1.0171e-137
GGGGTAACCATTCGAGTTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...
|||||||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGGGTAACCATTCGAGCTCATACAGATACACGGCTCCAAATCTCAAACGCTGAGAAAACAAAACCTCTCCTAATG...


# No relation to chimps did the engine search on NCBI Website Makes sense since this connected to a fruit fly.


```


# Open CV 1-3 

```




Code

Python 3
import numpy as np
import matplotlib.pyplot as plt 
%matplotlib inline 
import cv2
img = cv2.imread("German shepherd.jpeg")
type(img)
numpy.ndarray
img_wrong = cv2.imread('worong/path/doesnot?abcdegh/jpg')
type(img_wrong)
NoneType
plt.imshow(img)
<matplotlib.image.AxesImage at 0x7f652ed9abd0>

fix_img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
plt.imshow(fix_img)
<matplotlib.image.AxesImage at 0x7f652ecd6590>

img_gray = cv2.imread("German shepherd.jpeg", cv2.IMREAD_GRAYSCALE)
img_gray.shape
(176, 286)
plt.imshow(img_gray)
<matplotlib.image.AxesImage at 0x7f652d445610>

plt.imshow(img_gray, cmap = "gray")
<matplotlib.image.AxesImage at 0x7f652d3c1250>

fix_img.shape
(176, 286, 3)
new_image = cv2.resize(fix_img,(1000,400))
plt.imshow(new_image)
<matplotlib.image.AxesImage at 0x7f652c32d810>

new_image.shape
(400, 1000, 3)
w_ratio = 0.5
h_ratio = 0.5
​
new_image = cv2.resize(fix_img, (0,0), fix_img, w_ratio, h_ratio)
plt.imshow(new_image)
<matplotlib.image.AxesImage at 0x7f652c30d790>

new_image.shape
(88, 143, 3)
flip_img = cv2.flip(fix_img, 0)
plt.imshow(flip_img)
<matplotlib.image.AxesImage at 0x7f652c282c10>

flip_img2 = cv2.flip(fix_img, -1)
plt.imshow(flip_img2)
<matplotlib.image.AxesImage at 0x7f652c1f14d0>

type(fix_img)
numpy.ndarray
cv2.imwrite('German shepherd.jpeg', fix_img)
True

Code

Python 3
import cv2
import matplotlib.pyplot as plt 
%matplotlib inline 
img = cv2.imread("German shepherd.jpeg")
plt.imshow(img)
<matplotlib.image.AxesImage at 0x7fa69078b110>

img1 = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
plt.imshow(img1)
<matplotlib.image.AxesImage at 0x7fa690732f10>

img2 = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
plt.imshow(img2)
<matplotlib.image.AxesImage at 0x7fa68eead810>

img3 = cv2.cvtColor(img, cv2.COLOR_BGR2HLS)
plt.imshow(img3)
<matplotlib.image.AxesImage at 0x7fa68ee1eb90>

img1 = cv2.imread('Poodle.jpeg')
img2 = cv2.imread("German shepherd.jpeg")
plt.imshow(img1)
<matplotlib.image.AxesImage at 0x7fa68ed8ce90>

img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)
img2 = cv2.cvtColor(img2, cv2.COLOR_BGR2RGB)
plt.imshow(img1)
<matplotlib.image.AxesImage at 0x7fa68ed7bc10>

plt.imshow(img2)
<matplotlib.image.AxesImage at 0x7fa68eceab10>

img1 = cv2.resize(img1,(1200,1200))
img2 = cv2.resize(img2, (1200,1200))
alpha = 0.5 
beta = 0.5 
blended = cv2.addWeighted(img1, alpha, img2, beta, gamma=0)
plt.imshow(blended)
<matplotlib.image.AxesImage at 0x7fa68ec5ed50>

alpha = 0.8
beta = 0.2 
​
blended1 = cv2.addWeighted(img1, alpha, img2, beta, 0)
plt.imshow(blended1)
<matplotlib.image.AxesImage at 0x7fa68c34c410>

img1 = cv2.imread('Poodle.jpeg')
img2 = cv2.imread('German shepherd.jpeg')
​
img1 = cv2. cvtColor(img1, cv2.COLOR_BGR2RGB)
img2 = cv2. cvtColor(img2, cv2.COLOR_BGR2RGB)
​
img1 = cv2.resize(img1, (200,200))
img2 = cv2.resize(img2, (600,600))
large_img = img2
small_img = img1
​
x_offset = 0
y_offset = 0
​
x_end = x_offset + small_img.shape[1]
y_end = y_offset + small_img.shape[0]
​
small_img_resized = cv2.resize(small_img, (x_end - x_offset, y_end - y_offset))
​
large_img[y_offset:y_end, x_offset:x_end] = small_img_resized
​
plt.imshow(large_img)
plt.show()

large_img = img2
small_img = img1
​
x_offset = 0
y_offset = 0
​
x_end = x_offset + small_img.shape[1]
y_end = y_offset + small_img.shape[0]
​
large_img[y_offset:y_end, x_offset:x_end] = small_img
​
plt.imshow(large_img)
<matplotlib.image.AxesImage at 0x7fa68e9e9810>


# https//github.com/worklifesg/Python-for-Computer-Vision-with-OpenCV-and-Deep-Learning
import cv2 
import matplotlib.pyplot as plt
%matplotlib inline 
img = cv2.imread('rainbow.jpg')
plt.imshow(img)
<matplotlib.image.AxesImage at 0x7fe808ca51d0>
​
img = cv2.imread('rainbow.jpg', 0)
plt.imshow(img, cmap = 'gray')
<matplotlib.image.AxesImage at 0x7fe808c48790>

ret1, thresh1 = cv2.threshold(img, 127, 255, cv2.THRESH_BINARY)
ret1
127.0
plt.imshow(thresh1, cmap = "gray")
<matplotlib.image.AxesImage at 0x7fe8083b3910>

img2 = cv2.imread('rainbow.jpg', 0)
ret1, thresh1 = cv2.threshold(img2, 127, 255, cv2.THRESH_TRUNC)
plt.imshow(thresh1, cmap = "gray")
<matplotlib.image.AxesImage at 0x7fe808398a90>

img3 = cv2.imread('rainbow.jpg', 0)
ret1, thresh1 = cv2.threshold(img3, 127, 255, cv2.THRESH_TOZERO)
plt.imshow(thresh1, cmap = "gray")
<matplotlib.image.AxesImage at 0x7fe808305650>

img_r = cv2.imread('crossword.jpg', 0)
plt.imshow(img_r, cmap = 'gray')
<matplotlib.image.AxesImage at 0x7fe808262e50>

def show_pic(img):
    fig = plt.figure(figsize = (15,15))
    ax = fig.add_subplot(111)
    ax.imshow(img, cmap ='gray')
show_pic(img_r)

ret, th1 = cv2.threshold(img_r, 127, 255, cv2.THRESH_BINARY)
show_pic(th1)

ret, th1 = cv2.threshold(img_r, 200,255, cv2.THRESH_BINARY)
show_pic(th1)

th2 = cv2.adaptiveThreshold(img_r, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY, 11, 8)
show_pic(th2)

blended = cv2.addWeighted(src1 = th1, alpha =0.6,
                        src2 = th2, beta =0.4, gamma =0)
show_pic(blended)

th3 = cv2.adaptiveThreshold(img_r, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY, 11, 8)
​
blended = cv2.addWeighted(src1=th1,alpha=0.6,
                          src2=th2,beta= 0.4,gamma=0)
show_pic(blended)

```


# Aspect Detection / Corner/ Edge Detection 

```

import cv2
import numpy as np 
import matplotlib.pyplot as plt 
%matplotlib inline 
flat_chess = cv2.imread('Chessboard_green.png')
flat_chess = cv2.cvtColor(flat_chess, cv2.COLOR_BGR2RGB)
plt.imshow(flat_chess)
<matplotlib.image.AxesImage at 0x7ff3029db110>

gray_flat_chess = cv2.cvtColor(flat_chess, cv2.COLOR_BGR2GRAY)
plt.imshow(gray_flat_chess, cmap = "gray")
real_chess = cv2.imread("Chessboard.jpg")
real_chess = cv2.cvtColor(real_chess, cv2.COLOR_BGR2RGB)
plt.imshow(real_chess)
<matplotlib.image.AxesImage at 0x7ff3022ccd50>

gray_real_chess = cv2.cvtColor(real_chess, cv2.COLOR_BGR2GRAY)
plt.imshow(gray_real_chess, cmap = 'gray')
<matplotlib.image.AxesImage at 0x7ff301245ed0>

gray = np.float32(gray_flat_chess)
dst = cv2.cornerHarris(src = gray, blockSize = 2, ksize = 3, k = 0.04)
​
dst = cv2.dilate(dst, None)
flat_chess[dst>0.01*dst.max()] = [255,0,0]
​
plt.imshow(flat_chess)
<matplotlib.image.AxesImage at 0x7ff3011b5fd0>

gray = np.float32(gray_real_chess)
dst = cv2.cornerHarris(src = gray, blockSize =2, ksize=3,k=0.04)
dst = cv2.dilate(dst, None)
​
real_chess[dst>0.01*dst.max()] = [255, 0, 0]
​
plt.imshow(real_chess)
<matplotlib.image.AxesImage at 0x7ff301133d90>

#Shi_Tomasi Corner Detection 
​
corners = cv2.goodFeaturesToTrack(gray_flat_chess, 64, 0.01, 10)
corners = np.int0(corners)
​
for i in corners:
    x,y = i.ravel()
    cv2.circle(flat_chess, (x,y),3,(255,0,0), -1)
​
    plt.imshow(flat_chess)

corners = cv2.goodFeaturesToTrack(gray_real_chess, 100, 0.01, 10)
​
corners = np.int0(corners)
​
for i in corners:
    x, y = i.ravel()
    cv2.circle(real_chess, (x,y), 3, (0,255,0), -1)
    
plt.imshow(real_chess)
<matplotlib.image.AxesImage at 0x7ff3010056d0>




```



```



import cv2
import numpy as np 
import matplotlib.pyplot as plt 
%matplotlib inline 
img = cv2.imread("German shepherd.jpeg")
plt.imshow(img)
<matplotlib.image.AxesImage at 0x7fe78e3d0bd0>

edges = cv2.Canny(image = img, threshold1 = 127, threshold2 =127)
​
plt.imshow(edges)
<matplotlib.image.AxesImage at 0x7fe78ca92d10>

med_value = np.median(img)
med_value
177.0
lower = int(max(0, 0.7*med_value))
upper = int(min(255,1.3*med_value))
​
edges = cv2.Canny(img, threshold1 = lower, threshold2 =upper)
​
plt.imshow(edges)
<matplotlib.image.AxesImage at 0x7fe78c20bc50>

edges = cv2.Canny(image = img, threshold1 = lower, threshold2 = upper +100)
​
plt.imshow(edges)
<matplotlib.image.AxesImage at 0x7fe78c177bd0>

blurred_img = cv2.blur(img, ksize = (5,5))
​
edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper + 50)
​
plt.imshow(edges)
<matplotlib.image.AxesImage at 0x7fe78c055610>

blurred_img = cv2.blur(img, ksize = (5,5))
​
edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper + 100)
​
plt.imshow(edges)
<matplotlib.image.AxesImage at 0x7fe78c044990>

blurred_img = cv2.blur(img, ksize = (7,7))
​
edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper + 60)
​
plt.imshow(edges)
<matplotlib.image.AxesImage at 0x7fe786e6f210>

blurred_img = cv2.blur(img, ksize = (8,8))
​
edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper)
​
plt.imshow(edges)
<matplotlib.image.AxesImage at 0x7fe786dd5510>



```


# Feature Detection/ Feature Matching/ Object Detection 


```



Code

Python 3
import cv2 
import numpy as np 
import matplotlib.pyplot as plt 
%matplotlib inline 
def display(img, cmap ='gray'):
    fig =plt.figure(figsize = (12,10))
    ax = fig.add_subplot(111)
    ax.imshow(img,cmap = 'gray')
apple_jacks = cv2.imread("Apple_Jacks.jpg", 0)
display(apple_jacks)

cereals = cv2.imread('All_Cereal.jpg', 0)
display(cereals)

orb = cv2.ORB_create()
​
kp1,des1 = orb.detectAndCompute(apple_jacks, mask=None)
kp2,des2 = orb.detectAndCompute(cereals, mask=None)
bf = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck = True)
matches = bf.match(des1, des2)
matches = sorted(matches, key = lambda x:x.distance )
apple_jacks_matches = cv2.drawMatches(apple_jacks, kp1, cereals, kp2, matches[:25], None, flags =2)
display(apple_jacks_matches)

sift = cv2.SIFT_create()
kp1, des1 = sift.detectAndCompute(apple_jacks, None)
kp2, des2 = sift.detectAndCompute(cereals, None)
bf = cv2.BFMatcher()
matches = bf.knnMatch(des1, des2, k=2)
good = []
​
for match1, match2 in matches:
    if match1.distance < 0.75*match2.distance:
        good.append([match1])
        
print('Length of total matches:', len(matches))
print('Length of good matches:', len(good))
Length of total matches: 679
Length of good matches: 2
sift_matches = cv2.drawMatchesKnn(apple_jacks, kp1, cereals, kp2, good, None, flags =2)
display(sift_matches)

sift = cv2.SIFT_create()
​
kp1, des1 = sift.detectAndCompute(apple_jacks, None)
kp2, des2 = sift.detectAndCompute(cereals, None)
flann_index_KDtree = 0 
index_params = dict(algorithm=flann_index_KDtree, trees = 5)
search_params =dict(checks=50)
flann = cv2.FlannBasedMatcher(index_params, search_params)
​
matches = flann.knnMatch(des1, des2, k=2)
​
good = []
​
for match1, match2, in matches:
    if match1.distance < 0.75*match2.distance:
        good.append([match1])
flann_matches = cv2.drawMatchesKnn(apple_jacks, kp1, cereals, kp2, good, None, flags=0)
display(flann_matches)

sift = cv2.SIFT_create()
​
kp1, des1 = sift.detectAndCompute(apple_jacks, None)
kp2, des2 = sift.detectAndCompute(cereals, None)
flann_index_KDtree = 0
indx_params = dict(algorithms= flann_index_KDtree, trees = 5)
search_param = dict(checks = 50)
flann = cv2.FlannBasedMatcher(index_params, search_params)
​
matches = flann.knnMatch(des1, des2, k = 2)
matchesMask = [[0,0] for i in range(len(matches))]
for i, (match1, match2) in enumerate(matches):
    if match1.distance <0.75*match2.distance:
        matchesMask[i] = [1,0]
        
draw_params = dict(matchColor = (0,255,0),
                  singlePointColor = (255,0,0),
                  matchesMask = matchesMask,
                  flags =0)
flann_matches = cv2.drawMatchesKnn(apple_jacks, kp1, cereals, kp2, matches, None, **draw_params)
​
display(flann_matches)



```


```

Python 3
import cv2
import numpy as np 
import matplotlib.pyplot as plt 
%matplotlib inline 
full = cv2.imread('Training_Sunflower6.jpg')
full = cv2.cvtColor(full, cv2.COLOR_BGR2RGB)
plt.imshow(full)
<matplotlib.image.AxesImage at 0x7fd445220090>

test = cv2.imread('Sunflowers_helianthus_annuus.jpg')
test = cv2.cvtColor(test, cv2.COLOR_BGR2RGB)
plt.imshow(test)
<matplotlib.image.AxesImage at 0x7fd44542de50>

print('Test image shape:', full.shape)
print('Training image shape:', test.shape)
Test image shape: (225, 225, 3)
Training image shape: (1999, 3229, 3)
methods = ['cv2.TM_CCOEFF', 'cv2.TM_CCOEFF_NORMED', 'cv2.TM_CCORR', 'cv2.TM_CCORR_NORMED', 'cv2.TM_SQDIFF', 'cv2.TM_SQDIFF_NORMED']
for m in methods:
​
    test_copy = test.copy()
    method = eval(m)
    
    res = cv2.matchTemplate(test_copy, full, method)
    
    min_val, max_val, min_loc, max_loc = cv2.minMaxLoc(res)
    
    if method in [cv2.TM_SQDIFF, cv2.TM_SQDIFF_NORMED]:
        top_left = min_loc
    else:
        top_left = max_loc
        
    height, width, channels = full.shape 
    bottom_right = (top_left[0] + width, top_left[1] + height)
    
    cv2.rectangle(test_copy, top_left, bottom_right, (255,0,0),10)
    
    plt.subplot(121)
    plt.imshow(res)
    plt.title("Heatmap of template matching")
    plt.subplot(122)
    plt.imshow(test_copy)
    plt.title('Detection of template')
    
    plt.suptitle(m)
    
    plt.show()
    print('\n')
    print('\n')

```

# FIN 






























​



















    
   
