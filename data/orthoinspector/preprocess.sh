curl http://lbgi.fr/orthoinspector/dbquery_qfo/data/CSV/62.csv.gz -o 62.csv.gz
gunzip -c 62.csv.gz > 62.csv
cat 62.csv | awk -F '|' 'BEGIN {} {
    cele = $2;
    hsap = $4;
    desc = $5;
    printf("%s,%s|%s\n",cele,hsap,desc); }
    END {}' | grep ";Homo sapiens;" | sed "s/|.*//g;" > orthologs.csv
rm 62.csv.gz 62.csv