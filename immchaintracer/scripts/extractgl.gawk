BEGIN {
  # datafiles are tab-delimited
  FS = "\t"

  # pre-defined column label and index of sequence_id
  seqid_label = "sequence_id"
  seqid_idx = 0

  # if database filename or germline sequence_id is not specified
  # exit with error
  if( ARGV[1] == "" )
  {
    print "Database TSV file is not specified!" > "/dev/stderr"
    exit 1
  }
  if( seqid == "" )
  {
    print "Germline sequence id (seqid) is not specified!" > "/dev/stderr"
    exit 1
  }
}

BEGINFILE {
  dbfile = FILENAME

  # identify which columns contain sequence_id string
  getline < dbfile
  for(i=1; i<=NF; i++)
  {
    if( $i == seqid_label ) seqid_idx = i
  }
  # if not found, exit with error
  if( seqid_idx == 0  )
  {
    print "Header does not contain 'sequence_id' column!" > "/dev/stderr"
    exit 1
  }

  close( FILENAME )
}

# copy header
FNR == 1 {
  print $0
}
# read and process input database file
FNR > 1 {
  # if germline sequence is found -> extract all data
  if( $seqid_idx == seqid )   print $0
}