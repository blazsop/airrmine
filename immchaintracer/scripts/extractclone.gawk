BEGIN {
  # datafiles are tab-delimited
  FS = "\t"

  # pre-defined column labels and indices of sequence_id, clone_id fields
  clone_label = "clone_id"
  clone_idx = 0
  clone = ""
  seqid_label = "sequence_id"
  seqid_idx = 0

  # if database filename or germline sequence id is not specified
  # exit with error
  if( ARGV[1] == "" )
  {
    print "Database TSV file is not specified!" > "/dev/stderr"
    exit 1
  }
}

BEGINFILE {
  idx=1
  dbfile = FILENAME
  dbname = FILENAME

  # retrieve germline ID from filename
  sub( ".*/",  "", dbname )
  split( dbname, dbname_parts, "_" )
  glseqid = dbname_parts[5]

  # sanitiy check: does the filename contain the germline ID ?
  if( glseqid == 0 )
  {
    print "Filename does not contain the germline's ID!" > "/dev/stderr"
    exit 1
  }

  # identify which columns contain "sequence_id", "clone_id" data
  getline < dbfile
  for(i=1; i<=NF; i++)
  {
    if( $i == clone_label ) clone_idx = i
    if( $i == seqid_label ) seqid_idx = i
  }
  # if not found any, exit with error
  if( clone_idx == 0 || seqid_idx == 0 )
  {
    print "Header does not contain 'sequence_id' and/or 'clone_id' columns!" > "/dev/stderr"
    exit 1
  }

  # loop through file to find clone_id for sequence_id
  while( (getline < dbfile) == 1 )
    if( $seqid_idx == glseqid )
    {
      clone = $clone_idx
      break
    }
  close( FILENAME )
}

# read and process input database file skipping header
FNR > 1 {
  # if sequence belongs to the specified germline's clone but it is NOT the
  # germline itself -> extract it
  if( $clone_idx == clone )
    if( !($seqid_idx == glseqid) )
      print $0
}
