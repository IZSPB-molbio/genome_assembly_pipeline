import glob, os, sys, requests, base64, json
from snakemake import shell

sample        = snakemake.wildcards.sample
assembly_file = snakemake.input.assembly
outdir        = snakemake.params.outdir
output_json   = snakemake.output.rmlst_json
output_tab    = snakemake.output.rmlst_tab
threads       = snakemake.threads
log           = snakemake.log[0]

print("Hello!")

# with open(log, "w") as f:
# sys.stderr = sys.stdout = f
sys.stderr = sys.stdout = open(log, 'w')
uri = 'http://rest.pubmlst.org/db/pubmlst_rmlst_seqdef_kiosk/schemes/1/sequence'
#    with open(args.file, 'r') as x: 
#        fasta = x.read()
fasta = open(assembly_file, 'r').read()
payload = '{"base64":true,"details":true,"sequence":"' + base64.b64encode(fasta.encode()).decode() + '"}'
response = requests.post(uri, data=payload)
if response.status_code == requests.codes.ok:
    data = response.json()
    try: 
        data['taxon_prediction']
    except KeyError:
        # f.write("No match")
        print("No match")
        sys.exit(0)
    # This is for logging
    for match in data['taxon_prediction']:
            print("Rank: " + match['rank'] + "\n")
            print("Taxon: " + match['taxon'] + "\n")
            print("Support: " + str(match['support']) + "%\n")
            print("Taxonomy: " + match['taxonomy'] + "\n")
    # This is for tab output
    outhandle_output_tab = open(output_tab, 'w') 
    for match in data['taxon_prediction']:
        outhandle_output_tab.write("{sample}\t{rank}\t{support}\t{taxon}\t{taxonomy}\n".format(
                                    sample=sample,
                                    rank=match['rank'],
                                    support=match['support'],
                                    taxon=match['taxon'],
                                    taxonomy=match['taxonomy']
                                    ))
    outhandle_output_tab.close()
else:
    print(response.text)

with open(output_json, 'w') as outhandle_output_json:
    outhandle_output_json.write(json.dumps(data, indent=4))

#f.close()
