import csv

# Input data
data = """

CM024069.1	StringTie	transcript	9074088	9077408	1000	-	.	gene_id "MSTRG.839"; transcript_id "MSTRG.839.4"; 
CM024069.1	StringTie	exon	9074088	9075083	1000	-	.	gene_id "MSTRG.839"; transcript_id "MSTRG.839.4"; exon_number "1"; 
CM024069.1	StringTie	exon	9075275	9077408	1000	-	.	gene_id "MSTRG.839"; transcript_id "MSTRG.839.4"; exon_number "2"; 
CM024069.1	StringTie	transcript	9082951	9098707	1000	-	.	gene_id "VMungo0251G2940"; transcript_id "VMungo0251G2940.1"; ref_gene_id "VMungo0251G2940"; 
CM024069.1	StringTie	exon	9082951	9083016	1000	-	.	gene_id "VMungo0251G2940"; transcript_id "VMungo0251G2940.1"; exon_number "1"; ref_gene_id "VMungo0251G2940"; 
CM024069.1	StringTie	exon	9084171	9084474	1000	-	.	gene_id "VMungo0251G2940"; transcript_id "VMungo0251G2940.1"; exon_number "2"; ref_gene_id "VMungo0251G2940"; 
CM024069.1	StringTie	exon	9084758	9085004	1000	-	.	gene_id "VMungo0251G2940"; transcript_id "VMungo0251G2940.1"; exon_number "3"; ref_gene_id "VMungo0251G2940"; 
CM024069.1	StringTie	exon	9098535	9098707	1000	-	.	gene_id "VMungo0251G2940"; transcript_id "VMungo0251G2940.1"; exon_number "4"; ref_gene_id "VMungo0251G2940"; 

"""

# Parse lines
lines = data.strip().split('\n')

# Collect output
output = []
for line in lines:
    parts = line.strip().split('\t')
    if parts[2] == 'exon':
        seqid = parts[0]
        start = parts[3]
        end = parts[4]
        strand = parts[6]
        attr = parts[8]
        transcript_id = [field for field in attr.split(';') if 'transcript_id' in field][0].split('"')[1]
        output.append([seqid, transcript_id, start, end, strand])

# Save to TSV file
output_file = "D:\VM_lncRNA\VmLNC_Gene Structure\STK_coordinates.tsv"
with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(["ID", "source", "start", "end", "strand"])
    writer.writerows(output)

print(f"Output saved to: {output_file}")
