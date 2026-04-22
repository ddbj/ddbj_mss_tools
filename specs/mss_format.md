# Specification for MSS conversion script

## chromosomes.txt ファイルの形式 (egapx2mss)
5列で左から seq_id, type, seq_name, status, topology とする。
\#のある列はヘッダーで無視をする。
seq_id は FASTA ファイルのヘッダーと同じ。
type は chromosome, organelle, unplaced のいずれか。デフォルトでは unplaced。  
seq_name は 染色体ならその番号 (例: 1,2,3,4,I,II,V,X,Y)、organelle の場合には chloroplast/mitochondrionなど。chromosomeが1種類しかない場合およびunplacedの場合は空欄を許容する。unplacedの場合は空。デフォルトは空。
status は complete か partial のいずれか。デフォルトは partial
topology は linear か circular のいずれかでデフォルトは linear

chromosomes.txt に記載がない配列についてはすべてデフォルトとし、WGS としての登録扱いにする。
```
#seq_id	type	seq_name	status	topology
NC_024795.2_RagTag	chromosome	0	complete	
NW_022610937.1_RagTag	chromosome	1	partial	
ptg000034l	chromosome	1	partial	
ptg000046l	chromosome	3	complete	circular
ptg000120c	organelle	mitochondrion	complete	circular
ptg000123l	organelle	chloroplast	partial
```

## INFRASPECIFIC_NAME_MODIFIER
common.json の INFRASPECIFIC_NAME_MODIFIER で指定された文字列に等しい qualifier が個体を識別する情報としてFF_DEFINITION に使われる
たとえば INFRASPECIFIC_NAME_MODIFIER に strain が指定されていた場合、/strain qualifierの値が FF_DEFINITION  (例: "{organism} {source_identifier} DNA, chromosome {seq_name}, complete sequence") の {source_identifier} に使われる。


## source feature と FF_DEFITION
common.json の SOURCE に記載されたものをベースとし、配列のtype, statusに応じて変える。
source feature の ff_definition は次のように決められる  

1. WGS 登録かどうかのチェック
    すべての配列が unplaced の場合は WGS 登録とみなすので、以下のようにする  
    source feature に /submitter_seqid={seq_id}  
    ff_definition:  
    "{organism} {source_identifier} DNA, {seq_id}"    
    ここで seq_id は FASTA ヘッダーの ID　とする

2. 配列 type が chromosome の場合  
    以下で isoalate または strain は排他的に用いられる個体識別名。seq_name は 1, 2, 3 あるいは I, II、X, Y などの染色体名とする (chr1なら 1 だけを記載)
    2.1 status が complete の場合 ->   
        source feature に /chromosome qualifier を記載し、value に染色体番号を記載  
        ff_definition:  
        "{organism} {source_identifier} DNA, chromosome {seq_name}, complete sequence"
    2.2 status が complete 意外の場合 ->
        source feature に /chromosome qualifier を記載し、value に染色体番号を記載  
        ff_definition:  
        "{organism} {source_identifier} DNA, chromosome {seq_name}, unlocalized sequence {seq_id}"

    ただし、いずれの場合も seq_name が空 ("") の場合には /chromosome は不要で "chromosome {seq_name}" の部分は "chromosome" とする。
    complete の場合必ずしも T2T になっている必要はなく、全長の大部分を構築できていればギャップを含んでいても良いものとする。

3. 配列 type が organelle の場合
    3.1 status が complete の場合 ->   
        source feature に /organelle qualifier を記載し、value に染色体番号を記載  
        ff_definition:  
        "{organism} {source_identifier} DNA, {organelle_name}, complete sequence"
    3.2 status が complete 意外の場合 ->
        source feature に /chromosome qualifier を記載し、value に染色体番号を記載  
        ff_definition:  
        "{organism} {source_identifier} DNA, {organelle_name}, partial sequence"

4. 配列 type が unplaced の場合
    ff_definition: 
    "{organism} {source_identifier} DNA, unplaced sequence {seq_id}"
    とする。
