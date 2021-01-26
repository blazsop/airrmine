BEGIN {  RS = "\n"; FS = "\t" }
# print header once as the first line, then merge contents of all listed files
NR == 1 || FNR > 1 {  print  }
