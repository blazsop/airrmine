BEGIN {  RS = "\n"; FS = "," }
# print header once as the first line, then merge contents of all listed files
FNR > 1 || NR == 1 { print }
