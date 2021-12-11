for i in *faa; do
diamond blastp --query $i --db /scratch/1/arkadiy/databases/nr/nr.faa --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle --out ${i%.*}.nr.blast --max-target-seqs 1 --evalue 1E-6
done
