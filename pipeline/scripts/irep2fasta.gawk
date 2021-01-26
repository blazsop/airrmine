BEGIN {  FS=","  }
FNR > 1 {
print "> seq-" (FNR-1) "\n" $22
}
