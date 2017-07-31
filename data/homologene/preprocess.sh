gunzip -c homologene.data.gz > homologene.data
grep -E "^[0-9]*\t(6239|9606)" homologene.data > homologene.tsv
rm homologene.data
