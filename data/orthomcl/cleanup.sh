gunzip -c orthologs.txt.gz > orthologs.txt
grep "cele|" orthologs.txt | grep "hsap|" | sed '
	s/cele\|//g;
	s/	hsap\|/,/g;
	s/	.*//g;' > orthologs.csv

# grep "cele|" orthologs.txt | grep "hsap|" > OrthoMCL_cele-hsap.txt
rm orthologs.txt