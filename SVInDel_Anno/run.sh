python SVInDels_Feature_Annotation.py -g TAIR10_gene.gff -s ../final.gt.txt -m ID -c Parent -o test
python SVInDels_Features_Position.py TAIR10_gene.gff ../final.gt.txt Out_SVInDel
