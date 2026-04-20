batch_wgs_builder というツールを作成したい
目的は、複数のゲノムからWGS登録用ファイルを一括して作成する。
入力ファイルは common 部分を記載したjsonファイルと、各ゲノム固有のDBLINKやSOURCE情報を含んだsample_list.tsvファイル。
出力は各ゲノム配列に対する .ann ファイル と .fa ファイルのペア。
入力ファイル例は examples/batch_wgs_builder/common_example.json と examples/batch_wgs_builder/sample_list_WGS.tsv 
common_example.json ファイルには SUBMITTER や REFERENCE など全ゲノムファイルに共通のメタデータが含まれる
sample_listファイルはヘッダー1行目がannファイルでのFEATURE名に対応、二行目がqualifier名に対応する。１列目はゲノム塩基配列のFASTAファイルへのパスを示し、１行で1サンプルのメタデータ情報となっている。
メタデータはcommon jsonファイルの内容とsample_listの各行の内容を合わせてannファイルを作成する。
common jsonに全サンプルに共通するDBLINKやSOURCE featureの情報を含めて書くこともでき、その場合にはsample_listに記載されたもので上書きされる。
assembly_gap の記述方法や、annファイルの作成方法は基本的に mss_builder のものと同じとする。
submission_category に指定がなければ、KEYWORD は WGS と STANDARD_DRAFT を記載する。
annファイルのCOMMONにDATATYPEはWGSを指定。
submission_category にMAG-WGSを指定できる。その場合にはCOMMONのDIVISIONにENVを記載。
KEYWORDには ENV, WGS, STANDARD_DRAFT, Metagenome Assembled Genome, MAG の５つを記載。さらにmetagenome_sourceとenvironmental_sampleをsource featureに記載する。environmental_sampleはboolean qualifierなので値は空 (５列目が空)となる。

ff_definitionの推奨記載方は {organism} {source_modifier}, {entry}
souce feature にsubmitter_seqid qualifierを記載しその値は {entry} と同じとする。

WGSファイルの例: https://docs.google.com/spreadsheets/d/15gLGL5FMV8gRt46ezc2Gmb-R1NbYsIGMssB0MyHkcwE/edit?pli=1&gid=1134992157#gid=1134992157
MAGファイルの例: https://docs.google.com/spreadsheets/d/15gLGL5FMV8gRt46ezc2Gmb-R1NbYsIGMssB0MyHkcwE/edit?pli=1&gid=1453206143#gid=1453206143