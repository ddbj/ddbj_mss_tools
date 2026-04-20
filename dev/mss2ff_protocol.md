新たにmssフォーマット (.ann or .annt.tsvファイルと.fasta) からDDBJのFLAT FILEを生成するプログラム mss2ff を作成したい。
入力ファイルの例として　https://dfast.ddbj.nig.ac.jp/analysis/annotation/7d52e3dd-145e-44ac-b9fe-4d0d441098e9/submission/DDBJ.annt.tsv
と https://dfast.ddbj.nig.ac.jp/analysis/annotation/7d52e3dd-145e-44ac-b9fe-4d0d441098e9/submission/DDBJ.seq.fa が使用できます。
出力ファイルの例は https://dfast.ddbj.nig.ac.jp/analysis/download/7d52e3dd-145e-44ac-b9fe-4d0d441098e9/annotation.gbk です
exampleディレクトリにサブディレクトリ作成してダウンロードしてください。
ただし、この出力ファイルの例はgenbank形式ファイルなので厳密には正しくありません。
違い以下のようなものがあります
実際にDDBJから登録されている正しいDDBJ Flat Fileの例として　https://getentry.ddbj.nig.ac.jp/getentry/na/AP014680?filetype=text　を参考にしてください。これもexampleにダウンロードしてください
1. lineage情報として Bacteria; Bacillati; Bacillota; Bacilli; Lactobacillales; Lactobacillaceae; Paucilactobacillus. のような情報が入ります。lineageはNCBIのAPIなどを利用して作成してください。
2. DDBJでは登録者情報は1番目のREFERENCE情報に記載されます (通称Reference 1と呼ばれ、https://www.ddbj.nig.ac.jp/ddbj/flat-file.html#Reference1B　に詳細が記載されます)上記の例では2番目にDirect Submissionとして記載されているので正しくありません。
3. ##Genome-Assembly-Data-START## と ##Genome-Assembly-Data-END## に囲まれた部分はST_COMMENTと呼ばれるものです。これはCOMMENTの項目に書かれますが、ST_COMMENTの前に通常のコメントが挿入され、間に１行空行が挿入されます。ST_COMMENTのヘッダとフッタには tagset_id で指定されているものを使ってください。DDBJでは現在のところ Genome-Assembly-Data を Assembly-Data に対応しています
4. 塩基配列が記載されている ORIGIN という行の前に配列中のacgtの塩基数が "BASE COUNT       702912 a       437272 c       431103 g       706698 t" のような文字列で挿入されています
また、ファイル作成における諸注意として以下を留意してください
a. DIVISION情報はannファイルには含まれないのでコマンドオプションとして受け取れるように。ない場合には不明 (UNK) にしてください。
b. 登録日 (Reference中のSubmitted) およびファイル作成日 (1行目の日付) をオプションとして指定できるように。どちらも与えられない場合にはプログラム実行日を自動で指定してください
c. 塩基配列は小文字で入ります
d. ACCESSIONはFASTAファイルのIDを使用してください。VERSIONにはそれに .1 が加わったものとします。
e. COMMON featureやsource featureにメタ表記 (配列長を表すE、エントリ名を示す@@[entry]@@には実際の値を埋め込んでください。エントリ名はFASTAファイルのIDと同じものです。source featureのqualifierを示す @@[organism]@@などの値もsource featureに実際に記載されているものに変更してください。
f. source feature内のff_definitionはflat fileのDEFINITION行の内容を示します。メタ表記が含まれていた場合には実際の値を埋めて作成してください。
g. この例には含まれていませんが、annファイルにkeywordが含まれていた場合にはそれもflatfileに記載してください。たとえば WGS, STANDARD_DRAFT などのkeywordが指定される場合があります
h. 必要に応じてBiopythonのSeqRecord, SeqIOを使用してください。DDBJ形式はほぼGenBank形式と同じですが、GenBankの形式には BASE COUNT は存在しません。GenBank形式を出力した後に修正を加える方法で対応できるかもしれません。
i. pseudo qualifierや pseudogene qualifierが記載されていないCDS featureにはアミノ酸配列情報を translation qualifierとして含みます。annファイルには記載されていないのでflat fileに変換時に作成してください。ただしFASTAファイルが与えられていない場合には
j. アミノ酸配列への翻訳は https://raw.githubusercontent.com/nigyta/translate_with_exception/d3c382242f1372afb2b49b47a245ba8dcf548cf4/translate_with_transl_except.py を使ってください。transl_except が記載されている場合にも対応したものです。

REFERENCEについての微修正
https://www.ddbj.nig.ac.jp/ddbj/file-format.html#reference を参照
status が Unpublished の場合で year が記載されていなければ不要
status が In press の場合、journal と year を記載して
{journal} ({year}) In press
という形式
status が “Published” の場合は、Qualifier: journal, volume, start_pageが必要で end_page をオプショナルで指定可能
{journal} {volume}, (start_page)[-{end_page}] ({year})
という形式