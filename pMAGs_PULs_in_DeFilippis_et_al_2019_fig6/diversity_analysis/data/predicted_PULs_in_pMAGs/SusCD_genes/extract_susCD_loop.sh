for i in $(ls -1 ../*.puls.tsv); do

name=${i##*/}
filename=${name%*.puls.tsv}

grep -e "susC" $i | awk '{print $3}' > ${filename}_susC_ids.txt
grep -e "susD" $i | awk '{print $3}' > ${filename}_susD_ids.txt


perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' ${filename}_susC_ids.txt ../../proteins_MAGs_DeFillipis/${filename}.faa > ${filename}_susC.faa
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' ${filename}_susD_ids.txt ../../proteins_MAGs_DeFillipis/${filename}.faa > ${filename}_susD.faa;
done
